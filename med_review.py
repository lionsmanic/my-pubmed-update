import streamlit as st
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, timedelta
import time
import requests
import json

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 文獻智庫 v3.0", page_icon="🧬", layout="wide")

# --- Session State 初始化 ---
if 'email_content' not in st.session_state:
    st.session_state.email_content = ""
if 'analyzed_count' not in st.session_state:
    st.session_state.analyzed_count = 0
if 'run_analysis' not in st.session_state:
    st.session_state.run_analysis = False

# --- 核心函數：取得可用模型 ---
def get_available_models(api_key):
    url = f"https://generativelanguage.googleapis.com/v1beta/models?key={api_key}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            models = [m['name'].replace('models/', '') for m in data.get('models', []) 
                      if 'generateContent' in m.get('supportedGenerationMethods', [])]
            return models
        return []
    except: return []

# --- 側邊欄 ---
with st.sidebar:
    st.header("⚙️ 設定與模型")
    
    # 1. API Key
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    # 模型選擇
    selected_model_name = None
    if api_key:
        with st.spinner("偵測模型中..."):
            available_models = get_available_models(api_key)
        if available_models:
            default_ix = 0
            if 'gemini-1.5-flash' in available_models:
                default_ix = available_models.index('gemini-1.5-flash')
            elif 'gemini-pro' in available_models:
                default_ix = available_models.index('gemini-pro')
            selected_model_name = st.selectbox("✅ 選擇模型:", available_models, index=default_ix)
        else:
            st.error("❌ 無法取得模型 (請檢查 Key)")

    # 2. Email
    if 'EMAIL_ADDRESS' in st.secrets:
        user_email = st.secrets['EMAIL_ADDRESS']
    else:
        user_email = st.text_input("您的 Email", "lionsmanic@gmail.com")
    
    if 'EMAIL_PASSWORD' in st.secrets:
        email_password = st.secrets['EMAIL_PASSWORD']
    else:
        email_password = st.text_input("Gmail App Password", type="password")

    st.divider()
    
    # 3. 搜尋設定 (升級版)
    st.subheader("🔍 搜尋條件")
    
    # 預設選單
    KEYWORDS = {
        "🥚 婦癌 (Gyn Onc)": ["cervical cancer", "ovarian cancer", "endometrial cancer", "immunotherapy", "robotic surgery", "sarcoma"],
        "🌊 海扶刀 (HIFU)": ["HIFU", "high intensity focused ultrasound", "uterine leiomyoma", "adenomyosis", "fibroid"],
        "🧬 精準/其他": ["genetic test", "targeted therapy", "pembrolizumab", "bevacizumab"]
    }
    
    selected_cats = st.multiselect("📚 選擇主題類別", list(KEYWORDS.keys()), default=["🥚 婦癌 (Gyn Onc)"])
    base_keywords = []
    for cat in selected_cats:
        base_keywords.extend(KEYWORDS[cat])
    
    # 手動增加關鍵字
    custom_keywords_str = st.text_input("➕ 手動增加關鍵字 (用逗號隔開)", help="例如: TP53, recurrence, toxicity")
    if custom_keywords_str:
        custom_kws = [k.strip() for k in custom_keywords_str.split(",") if k.strip()]
        base_keywords.extend(custom_kws)
    
    final_keywords = st.multiselect("✅ 確認最終關鍵字", base_keywords, default=base_keywords)

    st.divider()
    
    # 期刊設定
    PRESET_JOURNALS = ["New England Journal of Medicine", "The Lancet", "The Lancet Oncology", "Journal of Clinical Oncology", "Gynecologic Oncology", "Journal of Gynecologic Oncology", "Nature Medicine"]
    use_journals = st.checkbox("限定期刊?", value=True)
    final_journals = []
    
    if use_journals:
        selected_journals = st.multiselect("選擇預設期刊", PRESET_JOURNALS, default=PRESET_JOURNALS)
        # 手動增加期刊
        custom_journals_str = st.text_input("➕ 手動增加期刊 (用逗號隔開)", help="例如: British Journal of Cancer")
        if custom_journals_str:
            custom_js = [j.strip() for j in custom_journals_str.split(",") if j.strip()]
            selected_journals.extend(custom_js)
        final_journals = selected_journals

    st.divider()

    # 4. 時間設定 (升級版：時間區段)
    st.subheader("📅 時間區段設定")
    date_mode = st.radio("時間模式", ["最近幾天", "指定過去區間"], index=0)
    
    date_range_query = ""
    
    if date_mode == "最近幾天":
        days_back = st.slider("搜尋過去幾天?", 1, 60, 7)
        # 使用 reldate 邏輯 (在函數中處理)
        date_params = {"reldate": days_back}
        display_date_info = f"過去 {days_back} 天"
    else:
        col1, col2 = st.columns(2)
        with col1:
            day_start = st.number_input("從幾天前開始?", min_value=1, value=60)
        with col2:
            day_end = st.number_input("到幾天前結束?", min_value=0, value=30)
        
        if day_start <= day_end:
            st.error("開始天數必須大於結束天數 (例如：從 60 天前 到 30 天前)")
            st.stop()
            
        # 計算日期字串 YYYY/MM/DD
        today = datetime.now()
        date_min = (today - timedelta(days=day_start)).strftime("%Y/%m/%d")
        date_max = (today - timedelta(days=day_end)).strftime("%Y/%m/%d")
        
        # PubMed 語法: "YYYY/MM/DD"[Date - Publication] : "YYYY/MM/DD"[Date - Publication]
        date_range_query = f' AND ("{date_min}"[Date - Publication] : "{date_max}"[Date - Publication])'
        date_params = {} # 這種模式下不用 reldate
        display_date_info = f"{date_min} ~ {date_max}"

    max_results = st.slider("篇數上限", 1, 10, 3)
    
    if st.button("🚀 開始搜尋與分析", type="primary", disabled=(not selected_model_name)):
        st.session_state.run_analysis = True
        st.session_state.email_content = ""
        st.session_state.analyzed_count = 0

# --- 核心功能 ---

def build_query(keywords, journals, date_str_query):
    if not keywords: return ""
    term_q = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    
    final = term_q
    if journals:
        journal_q = "(" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
        final = f"{term_q} AND {journal_q}"
    
    # 加上自訂的時間區間語法
    if date_str_query:
        final += date_str_query
        
    return final

def fetch_data(query, date_params, limit, email):
    Entrez.email = email
    try:
        # 如果是 reldate 模式，參數會放在 kwargs
        search_args = {
            "db": "pubmed",
            "term": query,
            "retmax": limit,
            "sort": "date"
        }
        if "reldate" in date_params:
            search_args["reldate"] = date_params["reldate"]
            
        h = Entrez.esearch(**search_args)
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
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    headers = {'Content-Type': 'application/json'}
    
    # 特殊指令：要求 AI 用 ||| 分隔「一句話簡介」與「詳細內容」
    prompt_text = f"""
    角色：資深婦科腫瘤醫師。
    任務：分析以下文獻。
    
    標題：{art['title']}
    摘要：{art['abstract']}
    
    【輸出格式要求 - 重要】：
    請輸出兩部分，中間用 "|||" 三個直槓符號嚴格區隔。
    
    第一部分：一句最精鍊的中文簡述 (One-liner)，告訴我這篇在做什麼，類似新聞標題。
    |||
    第二部分：詳細分析報告 (HTML 格式，不含 markdown)。內容須包含：
    1. 🧪 研究方法 (Methods): 簡述研究設計 (Retrospective? RCT? sample size?)
    2. 💡 發想緣起 (Rationale): 作者為何做這個？解決什麼臨床痛點？
    3. 📊 結果與數據 (Results): 重點 P 值、HR。
    4. 🏥 臨床運用與結論 (Conclusion): 婦癌醫師如何應用？
    
    HTML 樣式：使用 <div> <ul> <li> <b> 等標籤。
    """
    
    payload = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            result = response.json()
            try:
                full_text = result['candidates'][0]['content']['parts'][0]['text']
                # 分割字串
                parts = full_text.split("|||")
                if len(parts) >= 2:
                    summary = parts[0].strip()
                    detail_html = parts[1].strip().replace("```html", "").replace("```", "")
                    return summary, detail_html
                else:
                    return "分析完成", full_text # fallback
            except:
                return "解析錯誤", "<div style='color:red'>AI 回傳格式異常</div>"
        else:
            return "連線錯誤", f"<div style='color:red'>API 請求失敗: {response.text}</div>"
    except Exception as e:
        return "系統錯誤", f"<div style='color:red'>錯誤: {str(e)}</div>"

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
st.title("🧬 GynOnc 文獻智庫 v3.0")
st.caption("AI 驅動的婦科腫瘤精準文獻分析")

if st.session_state.run_analysis:
    if not api_key: st.warning("請輸入 API Key")
    elif not selected_model_name: st.warning("請選擇模型")
    else:
        with st.status("🔄 智能處理中...", expanded=True) as status:
            # 建構搜尋語法
            q = build_query(final_keywords, final_journals, date_range_query)
            st.write(f"📡 搜尋區間: `{display_date_info}`")
            # st.code(q) # debug用
            
            arts = fetch_data(q, date_params, max_results, user_email)
            
            if not arts:
                status.update(label="❌ 該時段無符合條件文章", state="error")
                st.session_state.run_analysis = False
            else:
                st.write(f"✅ 找到 {len(arts)} 篇，AI 正在深入閱讀...")
                st.session_state.email_content = ""
                cont = st.container()
                
                for i, art in enumerate(arts):
                    st.write(f"🤖 分析 #{i+1}: {art['title'][:30]}...")
                    
                    # 取得 簡述 和 詳細HTML
                    summary, detail_html = run_ai_direct_api(art, api_key, selected_model_name)
                    
                    with cont:
                        # 兩段式呈現
                        st.markdown("---")
                        # 第一段：標題 + 期刊 + 簡述
                        st.subheader(f"#{i+1} {art['title']}")
                        st.caption(f"📖 {art['journal']} | 🗓️ {display_date_info} | 🔗 [原文連結]({art['link']})")
                        st.info(f"📌 **精華速讀**: {summary}")
                        
                        # 第二段：展開看詳細
                        with st.expander("🩺 點擊查看：研究方法、發想緣起與詳細數據"):
                            st.markdown(detail_html, unsafe_allow_html=True)
                    
                    # Email 內容：標題 + 簡述 + 詳細
                    st.session_state.email_content += f"""
                    <div style="margin-bottom: 30px; border: 1px solid #ddd; padding: 15px; border-radius: 8px;">
                        <h3 style="color:#0056b3; margin-top:0;"><a href='{art['link']}'>{art['title']}</a></h3>
                        <p style="color:#666; font-size:0.9em;">{art['journal']}</p>
                        <div style="background:#eef6fc; padding:10px; border-radius:4px; margin-bottom:10px; color:#2c3e50; font-weight:bold;">
                            📌 {summary}
                        </div>
                        {detail_html}
                    </div>
                    """
                    time.sleep(1)
                
                st.session_state.analyzed_count = len(arts)
                status.update(label="🎉 分析完成！", state="complete")
                st.session_state.run_analysis = False

if st.session_state.analyzed_count > 0:
    st.divider()
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("📩 寄出彙整報告", type="primary"):
            if not email_password: st.error("需輸入 Gmail App Password")
            else:
                with st.spinner("寄信中..."):
                    ok, m = send_mail(user_email, email_password, st.session_state.email_content)
                    if ok: st.success(m)
                    else: st.error(m)
