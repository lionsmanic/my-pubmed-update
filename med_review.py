import streamlit as st
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
import time
import requests
import json

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 文獻智庫 (Auto-Detect)", page_icon="🧬", layout="wide")

# --- Session State 初始化 ---
if 'email_content' not in st.session_state:
    st.session_state.email_content = ""
if 'analyzed_count' not in st.session_state:
    st.session_state.analyzed_count = 0
if 'run_analysis' not in st.session_state:
    st.session_state.run_analysis = False

# --- 核心函數：取得可用模型清單 ---
def get_available_models(api_key):
    """直接詢問 Google 這把 Key 能用哪些模型"""
    url = f"https://generativelanguage.googleapis.com/v1beta/models?key={api_key}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            # 過濾出支援 generateContent 的模型
            models = [m['name'].replace('models/', '') for m in data.get('models', []) 
                      if 'generateContent' in m.get('supportedGenerationMethods', [])]
            return models
        else:
            return []
    except:
        return []

# --- 側邊欄 ---
with st.sidebar:
    st.header("⚙️ 設定與模型偵測")
    
    # 1. API Key
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    # --- 關鍵修改：模型選擇器 ---
    selected_model_name = None
    if api_key:
        with st.spinner("正在偵測可用模型..."):
            available_models = get_available_models(api_key)
        
        if available_models:
            # 優先尋找 flash 或 pro，否則選第一個
            default_ix = 0
            if 'gemini-1.5-flash' in available_models:
                default_ix = available_models.index('gemini-1.5-flash')
            elif 'gemini-pro' in available_models:
                default_ix = available_models.index('gemini-pro')
                
            selected_model_name = st.selectbox(
                "✅ 偵測到您的可用模型 (請選擇):", 
                available_models, 
                index=default_ix
            )
            st.caption(f"目前使用: {selected_model_name}")
        else:
            st.error("❌ 無法取得模型清單。請確認 API Key 是否正確，或是否已在 Google AI Studio 啟用 API。")
    # ---------------------------

    # 2. Email
    if 'EMAIL_ADDRESS' in st.secrets:
        user_email = st.secrets['EMAIL_ADDRESS']
    else:
        user_email = st.text_input("您的 Email", "lionsmanic@gmail.com")
        
    if 'EMAIL_PASSWORD' in st.secrets:
        email_password = st.secrets['EMAIL_PASSWORD']
    else:
        email_password = st.text_input("Gmail 應用程式密碼", type="password")

    st.divider()
    
    # 3. 搜尋設定
    KEYWORDS = {
        "🥚 婦癌 (Gyn Onc)": [
            "cervical cancer", "ovarian cancer", "endometrial cancer", 
            "immunotherapy", "robotic surgery", "sarcoma", 
            "gynecologic neoplasms"
        ],
        "🌊 海扶刀 (HIFU)": [
            "HIFU", "high intensity focused ultrasound", 
            "uterine leiomyoma", "adenomyosis", "fibroid"
        ],
        "🧬 其他/精準醫療": [
            "genetic test", "targeted therapy"
        ]
    }
    
    JOURNALS = [
        "New England Journal of Medicine", "Nature", "Science", "Cell", 
        "The Lancet", "The Lancet Oncology", "Nature Communications", 
        "Journal of Clinical Oncology", "JAMA", 
        "Gynecologic Oncology", "Journal of Gynecologic Oncology"
    ]

    st.subheader("🔍 搜尋參數")
    selected_categories = st.multiselect("選擇類別", list(KEYWORDS.keys()), default=["🥚 婦癌 (Gyn Onc)"])
    active_keywords = []
    for cat in selected_categories:
        active_keywords.extend(KEYWORDS[cat])
    final_keywords = st.multiselect("微調關鍵字", active_keywords, default=active_keywords)

    use_specific_journals = st.checkbox("限定權威期刊?", value=True)
    if use_specific_journals:
        selected_journals = st.multiselect("選擇期刊", JOURNALS, default=JOURNALS)
    
    days_back = st.slider("搜尋過去幾天?", 1, 60, 7)
    max_results = st.slider("篇數上限", 1, 10, 3)
    
    # 只有當選到了模型才允許執行
    if st.button("🚀 開始搜尋與分析", type="primary", disabled=(not selected_model_name)):
        st.session_state.run_analysis = True
        st.session_state.email_content = ""
        st.session_state.analyzed_count = 0

# --- 核心功能 ---

def build_query(keywords, journals):
    if not keywords: return ""
    term_q = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    if journals:
        journal_q = "(" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
        return f"{term_q} AND {journal_q}"
    return term_q

def fetch_data(query, days, limit, email):
    Entrez.email = email
    try:
        h = Entrez.esearch(db="pubmed", term=query, reldate=days, retmax=limit, sort="date")
        r = Entrez.read(h)
        ids = r["IdList"]
        if not ids: return []
        
        h = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        arts = Entrez.read(h)
        parsed = []
        for art in arts['PubmedArticle']:
            try:
                cit = art['MedlineCitation']
                ti = cit['Article']['ArticleTitle']
                jo = cit['Article']['Journal']['Title']
                ab = " ".join([str(x) for x in cit['Article']['Abstract']['AbstractText']]) if 'Abstract' in cit['Article'] else "No Abstract"
                ids = art['PubmedData']['ArticleIdList']
                doi = next((i for i in ids if i.attributes['IdType']=='doi'), None)
                link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{ids[0]}/"
                parsed.append({"title":ti, "journal":jo, "abstract":ab, "link":link})
            except: continue
        return parsed
    except Exception as e:
        st.error(f"PubMed Error: {e}"); return []

def run_ai_direct_api(art, key, model_name):
    """
    使用使用者選擇的模型名稱進行分析
    """
    # 這裡使用選擇的模型名稱
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    
    headers = {'Content-Type': 'application/json'}
    
    prompt_text = f"""
    角色：婦科腫瘤專家。請將以下摘要轉成繁體中文臨床重點 (HTML)。
    
    標題：{art['title']}
    摘要：{art['abstract']}
    
    輸出 HTML (不含markdown, 僅內容):
    <div style="background:#f9f9f9; padding:15px; border-left:4px solid #007bff; margin-bottom:10px;">
        <h4 style="color:#0056b3; margin-top:0;">📝 重點摘要</h4>
        <ul>
            <li><b>背景</b>: ...</li>
            <li><b>結果</b>: (含數據)...</li>
            <li><b>結論</b>: ...</li>
        </ul>
        <h4 style="color:#d35400;">💡 臨床洞察</h4>
        <ul>
            <li><b>發想緣起</b>: ...</li>
            <li><b>臨床運用</b>: ...</li>
            <li><b>未來機會</b>: ...</li>
        </ul>
    </div>
    """
    
    payload = {
        "contents": [{
            "parts": [{"text": prompt_text}]
        }]
    }
    
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        
        if response.status_code == 200:
            result = response.json()
            try:
                return result['candidates'][0]['content']['parts'][0]['text']
            except:
                return f"<div style='color:red'>❌ AI 回傳格式異常</div>"
        else:
            return f"<div style='color:red'>❌ API 請求失敗 (Code {response.status_code}): {response.text}</div>"
            
    except Exception as e:
        return f"<div style='color:red'>❌ 連線錯誤: {str(e)}</div>"

def send_mail(to, pwd, html):
    msg = MIMEMultipart()
    msg['From'] = to
    msg['To'] = to
    msg['Subject'] = f"GynOnc Report {datetime.now().strftime('%Y-%m-%d')}"
    
    full_html = f"<html><body style='font-family:Arial;'>{html}</body></html>"
    msg.attach(MIMEText(full_html, 'html'))
    try:
        s = smtplib.SMTP('smtp.gmail.com', 587)
        s.starttls()
        s.login(to, pwd)
        s.send_message(msg); s.quit()
        return True, "已寄出"
    except Exception as e: return False, str(e)

# --- 主程式 ---
st.title("🧬 GynOnc 文獻智庫 (Auto-Select)")

if st.session_state.run_analysis:
    if not api_key: st.warning("請輸入 API Key")
    elif not selected_model_name: st.warning("⚠️ 尚未選擇模型，請檢查側邊欄")
    elif not final_keywords: st.warning("請選擇關鍵字")
    else:
        with st.status("🔄 處理中...", expanded=True) as status:
            q = build_query(final_keywords, selected_journals if use_specific_journals else None)
            st.write(f"搜尋: `{q[:60]}...`")
            arts = fetch_data(q, days_back, max_results, user_email)
            
            if not arts:
                status.update(label="❌ 無新文章", state="error")
                st.session_state.run_analysis = False
            else:
                st.write(f"✅ 找到 {len(arts)} 篇")
                st.session_state.email_content = ""
                cont = st.container()
                
                for i, art in enumerate(arts):
                    st.write(f"🤖 分析 #{i+1}...")
                    
                    # 傳入使用者選擇的模型
                    ai_html = run_ai_direct_api(art, api_key, selected_model_name)
                    
                    with cont:
                        st.subheader(f"{i+1}. {art['title']}")
                        st.caption(f"{art['journal']} | [Link]({art['link']})")
                        st.markdown(ai_html, unsafe_allow_html=True)
                        st.divider()
                    
                    st.session_state.email_content += f"<h3><a href='{art['link']}'>{art['title']}</a></h3><p>{art['journal']}</p>{ai_html}<hr>"
                    time.sleep(1)
                
                st.session_state.analyzed_count = len(arts)
                status.update(label="🎉 完成", state="complete")
                st.session_state.run_analysis = False

if st.session_state.analyzed_count > 0:
    if st.button("📩 寄出報告", type="primary"):
        if not email_password: st.error("需輸入 Gmail App Password")
        else:
            ok, m = send_mail(user_email, email_password, st.session_state.email_content)
            if ok: st.success(m)
            else: st.error(m)
