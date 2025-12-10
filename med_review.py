import streamlit as st
import google.generativeai as genai
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
import time

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
    "New England Journal of Medicine", 
    "Nature", 
    "Science", 
    "Cell", 
    "The Lancet", 
    "The Lancet Oncology", 
    "Nature Communications", 
    "Journal of Clinical Oncology", 
    "JAMA", 
    "Gynecologic Oncology", 
    "Journal of Gynecologic Oncology"
]

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 文獻智庫", page_icon="🧬", layout="wide")

# --- Session State 初始化 ---
if 'email_content' not in st.session_state:
    st.session_state.email_content = ""
if 'analyzed_count' not in st.session_state:
    st.session_state.analyzed_count = 0
if 'run_analysis' not in st.session_state:
    st.session_state.run_analysis = False

# --- 側邊欄：設定控制台 ---
with st.sidebar:
    st.header("⚙️ 設定控制台")
    
    # 1. API Key 設定 (優先讀取 Secrets)
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 Gemini API Key 已載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    # 2. Email 設定
    if 'EMAIL_ADDRESS' in st.secrets:
        user_email = st.secrets['EMAIL_ADDRESS']
    else:
        user_email = st.text_input("您的 Email", "lionsmanic@gmail.com")
        
    if 'EMAIL_PASSWORD' in st.secrets:
        email_password = st.secrets['EMAIL_PASSWORD']
        st.success("🔑 Gmail 密碼已載入")
    else:
        email_password = st.text_input("Gmail 應用程式密碼", type="password", help="若只需瀏覽不需寄信可不填")

    st.divider()
    
    # 3. 搜尋條件
    st.subheader("🔍 1. 選擇搜尋主題")
    selected_categories = st.multiselect("選擇類別", list(KEYWORDS.keys()), default=["🥚 婦癌 (Gyn Onc)"])
    
    active_keywords = []
    for cat in selected_categories:
        active_keywords.extend(KEYWORDS[cat])
    
    final_keywords = st.multiselect("微調搜尋關鍵字", active_keywords, default=active_keywords)

    st.subheader("📚 2. 期刊篩選")
    use_specific_journals = st.checkbox("限定於指定權威期刊?", value=True)
    if use_specific_journals:
        selected_journals = st.multiselect("選擇期刊", JOURNALS, default=JOURNALS)
    
    st.subheader("📅 3. 其他條件")
    days_back = st.slider("搜尋過去幾天?", 1, 60, 7)
    max_results = st.slider("分析篇數上限", 1, 10, 3)
    
    # 啟動按鈕
    if st.button("🚀 開始搜尋與分析", type="primary"):
        st.session_state.run_analysis = True
        # 重置之前的結果
        st.session_state.email_content = ""
        st.session_state.analyzed_count = 0

# --- 核心功能函數 ---

def build_pubmed_query(keywords, journals, days_back):
    if not keywords: return ""
    term_query = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    
    if journals:
        journal_query = "(" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
        final_query = f"{term_query} AND {journal_query}"
    else:
        final_query = term_query
    return final_query

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

def gemini_analyze(article, key):
    # 設定指定的最新模型
    target_model = 'gemini-1.5-flash'
    
    try:
        genai.configure(api_key=key)
        model = genai.GenerativeModel(target_model)
        
        prompt = f"""
        角色：你是婦科腫瘤學的資深臨床醫師與研究員。
        任務：閱讀以下文獻摘要，並轉化為繁體中文的臨床簡報。
        格式：請直接輸出 HTML 代碼 (不要包含 ```html 標籤)，以便嵌入網頁與郵件。
        
        文獻標題：{article['title']}
        期刊：{article['journal']}
        摘要：{article['abstract']}
        
        請依照以下結構撰寫 HTML (請使用 <div> 區塊包覆)：
        <div style="font-family: sans-serif; padding: 10px; border-left: 4px solid #3498db; background-color: #f8f9fa; border-radius: 4px;">
            <h4 style="color: #2c3e50; margin-top: 0;">1. 📝 重點摘要</h4>
            <ul style="color: #333;">
                <li><b>背景/目的</b>: ...</li>
                <li><b>結果 (數據)</b>: (請務必保留 P值、HR、OR 等重要統計數據)...</li>
                <li><b>結論</b>: ...</li>
            </ul>
            <h4 style="color: #e67e22; margin-bottom: 5px;">2. 💡 臨床洞察與發想</h4>
            <ul style="color: #333;">
                <li><b>發想來源</b>: (這篇文章是基於什麼臨床痛點或未解之謎？)</li>
                <li><b>臨床可行運用</b>: (對婦科腫瘤醫師而言，這改變了什麼處置流程？)</li>
                <li><b>未來研究機會</b>: (我們是否能模仿此研究？或有哪些延伸題目適合繼續發展？)</li>
            </ul>
        </div>
        """
        
        response = model.generate_content(prompt)
        return response.text

    except Exception as e:
        # 詳細錯誤處理
        error_msg = str(e)
        if "404" in error_msg:
            return f"""
            <div style="color: red; border: 1px solid red; padding: 10px; background: #fff0f0;">
                <b>❌ 模型載入失敗 ({target_model})</b><br>
                錯誤代碼 404 表示您的 Python 環境套件版本過舊。<br>
                請務必確認 GitHub 上的 <code>requirements.txt</code> 內容為 <code>google-generativeai>=0.8.3</code>，並重啟 Streamlit App。
            </div>
            """
        else:
            return f"<div style='color:red'>❌ AI 分析發生未預期錯誤: {error_msg}</div>"

def send_email_via_gmail(to_email, password, html_content):
    msg = MIMEMultipart()
    msg['From'] = to_email
    msg['To'] = to_email
    msg['Subject'] = f"GynOnc 文獻彙整報告 ({datetime.now().strftime('%Y-%m-%d')})"
    
    full_html = f"""
    <html>
    <body style="font-family: Arial, sans-serif; line-height: 1.6; color: #333;">
        <h2 style="color: #2c3e50;">🧬 婦科腫瘤文獻智慧報告</h2>
        <p>生成時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        <hr style="border: 1px solid #eee;">
        {html_content}
        <br>
        <p style="font-size: 0.8em; color: #999;">本郵件由 Streamlit AI 助手自動生成。</p>
    </body>
    </html>
    """
    msg.attach(MIMEText(full_html, 'html'))
    
    try:
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login(to_email, password)
        server.send_message(msg)
        server.quit()
        return True, "寄送成功！"
    except Exception as e:
        return False, f"寄送失敗: {e}"

# --- 主程式邏輯 ---

st.title("🧬 GynOnc 婦癌文獻智庫 (AI Assistant)")

# 執行分析
if st.session_state.run_analysis:
    if not api_key:
        st.warning("⚠️ 請在側邊欄輸入 Gemini API Key 才能開始。")
    elif not final_keywords:
        st.warning("⚠️ 請至少選擇一個搜尋關鍵字。")
    else:
        # 1. 搜尋
        with st.status("🔄 正在搜尋 PubMed...", expanded=True) as status:
            q = build_pubmed_query(final_keywords, selected_journals if use_specific_journals else None, days_back)
            st.write(f"搜尋語法: `{q[:100]}...`") 
            
            articles = fetch_pubmed(q, days_back, max_results, user_email)
            
            if not articles:
                status.update(label="❌ 最近沒有符合條件的新文章。", state="error")
                st.session_state.run_analysis = False 
            else:
                st.write(f"✅ 找到 {len(articles)} 篇，AI 正在逐篇閱讀分析...")
                
                # 清空並準備 Email 內容容器
                st.session_state.email_content = ""
                
                results_container = st.container()
                
                for i, art in enumerate(articles):
                    st.write(f"🤖 分析第 {i+1} 篇: {art['title'][:30]}...")
                    analysis_html = gemini_analyze(art, api_key)
                    
                    # 畫面顯示
                    with results_container:
                        st.markdown("---")
                        st.subheader(f"#{i+1} {art['title']}")
                        st.caption(f"📖 {art['journal']} | 🗓️ {days_back}天內 | 🔗 [原文連結]({art['link']})")
                        st.markdown(analysis_html, unsafe_allow_html=True)
                    
                    # Email 內容堆疊
                    st.session_state.email_content += f"""
                    <div style="margin-bottom: 30px; padding: 15px; background-color: #ffffff; border: 1px solid #ddd; border-radius: 5px;">
                        <h3 style="margin-top: 0; color: #1a5276;"><a href="{art['link']}" style="text-decoration: none; color: #1a5276;">{art['title']}</a></h3>
                        <p style="font-size: 0.9em; color: #666;">📖 {art['journal']}</p>
                        {analysis_html}
                    </div>
                    """
                    time.sleep(1) # 避免 API 呼叫過快
                
                st.session_state.analyzed_count = len(articles)
                status.update(label="🎉 分析完成！請查看下方結果或寄出郵件。", state="complete")
                st.session_state.run_analysis = False

# 寄信按鈕區
if st.session_state.analyzed_count > 0:
    st.divider()
    st.markdown("### 📧 彙整與分享")
    st.info("如果您滿意上方的分析結果，點擊下方按鈕將其寄到您的信箱。")
    
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("📩 立即寄出彙整報告", type="primary"):
            if not email_password:
                st.error("❌ 尚未設定 Gmail 應用程式密碼，無法寄信。請在側邊欄輸入。")
            else:
                with st.spinner("正在寄信中..."):
                    success, msg = send_email_via_gmail(user_email, email_password, st.session_state.email_content)
                    if success:
                        st.success(f"✅ {msg} 請檢查收件匣！")
                    else:
                        st.error(msg)
