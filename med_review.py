import streamlit as st
import google.generativeai as genai
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
import time
import importlib.metadata # 用來檢查套件版本

# --- 1. 預定義的專業關鍵字與期刊清單 ---
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

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 文獻智庫 (除錯版)", page_icon="🛠️", layout="wide")

# --- Session State 初始化 ---
if 'email_content' not in st.session_state:
    st.session_state.email_content = ""
if 'analyzed_count' not in st.session_state:
    st.session_state.analyzed_count = 0
if 'run_analysis' not in st.session_state:
    st.session_state.run_analysis = False

# --- 側邊欄：設定與診斷 ---
with st.sidebar:
    st.header("🛠️ 系統診斷控制台")
    
    # --- 顯示套件版本 (關鍵檢查點) ---
    try:
        pkg_version = importlib.metadata.version('google-generativeai')
        st.info(f"📦 套件版本: {pkg_version}\n(應 >= 0.7.0)")
    except:
        st.error("❌ 無法偵測套件版本，安裝可能不完整")

    st.divider()

    # 1. API Key 設定
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已從 Secrets 載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    # 2. Email 設定
    if 'EMAIL_ADDRESS' in st.secrets:
        user_email = st.secrets['EMAIL_ADDRESS']
    else:
        user_email = st.text_input("您的 Email", "lionsmanic@gmail.com")
        
    if 'EMAIL_PASSWORD' in st.secrets:
        email_password = st.secrets['EMAIL_PASSWORD']
    else:
        email_password = st.text_input("Gmail App Password", type="password")

    st.divider()
    
    # 3. 測試按鈕 (新增功能：不用搜尋，直接測連線)
    if st.button("🔌 測試 AI 連線 (不搜尋)"):
        if not api_key:
            st.error("請先輸入 API Key")
        else:
            try:
                genai.configure(api_key=api_key)
                model = genai.GenerativeModel('gemini-1.5-flash')
                res = model.generate_content("Hello, this is a test.")
                st.success(f"✅ 連線成功！AI 回應: {res.text}")
            except Exception as e:
                st.error(f"❌ 連線失敗: {e}")

    st.subheader("🔍 搜尋設定")
    selected_categories = st.multiselect("選擇類別", list(KEYWORDS.keys()), default=["🥚 婦癌 (Gyn Onc)"])
    active_keywords = []
    for cat in selected_categories:
        active_keywords.extend(KEYWORDS[cat])
    final_keywords = st.multiselect("微調關鍵字", active_keywords, default=active_keywords)

    use_specific_journals = st.checkbox("限定權威期刊?", value=True)
    if use_specific_journals:
        selected_journals = st.multiselect("選擇期刊", JOURNALS, default=JOURNALS)
    
    days_back = st.slider("過去幾天?", 1, 60, 7)
    max_results = st.slider("篇數上限", 1, 10, 3)
    
    if st.button("🚀 開始搜尋與分析", type="primary"):
        st.session_state.run_analysis = True
        st.session_state.email_content = ""
        st.session_state.analyzed_count = 0

# --- 核心功能函數 ---

def build_pubmed_query(keywords, journals, days_back):
    if not keywords: return ""
    term_query = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    if journals:
        journal_query = "(" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
        return f"{term_query} AND {journal_query}"
    return term_query

def fetch_pubmed(query, days, max_res, email):
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="pubmed", term=query, reldate=days, retmax=max_res, sort="date")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        articles = Entrez.read(handle)
        
        parsed = []
        for art in articles['PubmedArticle']:
            try:
                cit = art['MedlineCitation']
                title = cit['Article']['ArticleTitle']
                journal = cit['Article']['Journal']['Title']
                abstract = " ".join([str(x) for x in cit['Article']['Abstract']['AbstractText']]) if 'Abstract' in cit['Article'] else "No Abstract"
                ids = art['PubmedData']['ArticleIdList']
                doi = next((item for item in ids if item.attributes['IdType'] == 'doi'), None)
                link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{ids[0]}/"
                parsed.append({"title": title, "journal": journal, "abstract": abstract, "link": link})
            except: continue
        return parsed
    except Exception as e:
        st.error(f"PubMed 連線錯誤: {e}")
        return []

def gemini_analyze_debug(article, key):
    # --- 這是除錯版本的分析函數 ---
    genai.configure(api_key=key)
    
    # 嘗試列表
    models = ['gemini-1.5-flash', 'gemini-1.5-flash-latest', 'gemini-pro']
    last_error = ""
    
    for model_name in models:
        try:
            model = genai.GenerativeModel(model_name)
            
            prompt = f"""
            角色：婦科腫瘤專家。
            任務：將以下摘要轉為繁體中文臨床重點 (HTML格式)。
            
            標題：{article['title']}
            摘要：{article['abstract']}
            
            請輸出 HTML (不含markdown標記)：
            <div style="font-family: sans-serif;">
                <h4 style="color: #2e86c1;">1. 📝 重點摘要 ({model_name})</h4>
                <ul>
                    <li><b>背景</b>: ...</li>
                    <li><b>結果</b>: (含 P值/HR)...</li>
                    <li><b>結論</b>: ...</li>
                </ul>
                <h4 style="color: #d35400;">2. 💡 臨床洞察</h4>
                <ul>
                    <li><b>發想來源</b>: ...</li>
                    <li><b>臨床運用</b>: ...</li>
                    <li><b>未來機會</b>: ...</li>
                </ul>
            </div>
            """
            response = model.generate_content(prompt)
            return response.text # 成功就回傳
        except Exception as e:
            last_error = str(e)
            continue # 失敗就試下一個

    # 如果全部失敗，回傳詳細錯誤
    error_html = f"""
    <div style="background: #ffe6e6; padding: 10px; border: 1px solid red; color: red;">
        <h3>⚠️ AI 分析失敗</h3>
        <p><b>錯誤原因:</b> {last_error}</p>
        <p><b>建議:</b> 請檢查 API Key 是否正確，或 requirements.txt 版本是否 >=0.7.0</p>
    </div>
    """
    return error_html

def send_email_via_gmail(to_email, password, html_content):
    msg = MIMEMultipart()
    msg['From'] = to_email
    msg['To'] = to_email
    msg['Subject'] = f"GynOnc 文獻報告 (Debug Mode) - {datetime.now().strftime('%Y-%m-%d')}"
    
    full_html = f"<html><body>{html_content}</body></html>"
    msg.attach(MIMEText(full_html, 'html'))
    
    try:
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login(to_email, password)
        server.send_message(msg)
        server.quit()
        return True, "寄送成功"
    except Exception as e:
        return False, f"寄送失敗: {e}"

# --- 主程式邏輯 ---

st.title("🧬 GynOnc 婦癌文獻智庫 (除錯版)")

if st.session_state.run_analysis:
    if not api_key:
        st.warning("⚠️ 請輸入 Gemini API Key")
    elif not final_keywords:
        st.warning("⚠️ 請選擇關鍵字")
    else:
        with st.status("🔄 執行中...", expanded=True) as status:
            q = build_pubmed_query(final_keywords, selected_journals if use_specific_journals else None, days_back)
            st.write(f"搜尋語法: `{q[:50]}...`")
            
            articles = fetch_pubmed(q, days_back, max_results, user_email)
            
            if not articles:
                status.update(label="❌ 找不到文章", state="error")
                st.session_state.run_analysis = False
            else:
                st.write(f"✅ 找到 {len(articles)} 篇")
                st.session_state.email_content = ""
                results_container = st.container()
                
                for i, art in enumerate(articles):
                    st.write(f"🤖 分析第 {i+1} 篇: {art['title'][:30]}...")
                    
                    # 呼叫除錯版函數
                    analysis_html = gemini_analyze_debug(art, api_key)
                    
                    with results_container:
                        st.markdown("---")
                        st.subheader(f"#{i+1} {art['title']}")
                        st.caption(f"📖 {art['journal']} | 🔗 [連結]({art['link']})")
                        st.markdown(analysis_html, unsafe_allow_html=True)
                    
                    st.session_state.email_content += f"""
                    <div style="margin-bottom:20px; padding:15px; background:#f9f9f9;">
                        <h3><a href="{art['link']}">{art['title']}</a></h3>
                        <p>{art['journal']}</p>
                        {analysis_html}
                    </div>
                    """
                    time.sleep(1)
                
                st.session_state.analyzed_count = len(articles)
                status.update(label="🎉 分析完成", state="complete")
                st.session_state.run_analysis = False

if st.session_state.analyzed_count > 0:
    st.divider()
    if st.button("📩 寄出測試郵件", type="primary"):
        if not email_password:
            st.error("請輸入 Gmail 應用程式密碼")
        else:
            success, msg = send_email_via_gmail(user_email, email_password, st.session_state.email_content)
            if success: st.success(msg)
            else: st.error(msg)
