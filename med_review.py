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
st.set_page_config(page_title="GynOnc 文獻快篩系統 v4.0", page_icon="⚡", layout="wide")

# --- Session State 初始化 (關鍵：用於記住搜尋結果與分析緩存) ---
if 'articles_data' not in st.session_state:
    st.session_state.articles_data = []  # 存搜尋到的所有文章簡介
if 'analysis_cache' not in st.session_state:
    st.session_state.analysis_cache = {} # 存已經分析過的詳細內容 {doi: html}
if 'email_queue' not in st.session_state:
    st.session_state.email_queue = []    # 準備寄出的清單

# --- 核心函數：模型偵測 ---
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
    st.header("⚡ 設定與快篩")
    
    # 1. API Key
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    # 模型選擇
    selected_model_name = None
    if api_key:
        available_models = get_available_models(api_key)
        if available_models:
            default_ix = 0
            if 'gemini-1.5-flash' in available_models: default_ix = available_models.index('gemini-1.5-flash')
            elif 'gemini-pro' in available_models: default_ix = available_models.index('gemini-pro')
            selected_model_name = st.selectbox("✅ AI 模型:", available_models, index=default_ix)

    # 2. Email
    if 'EMAIL_ADDRESS' in st.secrets:
        user_email = st.secrets['EMAIL_ADDRESS']
    else:
        user_email = st.text_input("Email", "lionsmanic@gmail.com")
    
    if 'EMAIL_PASSWORD' in st.secrets:
        email_password = st.secrets['EMAIL_PASSWORD']
    else:
        email_password = st.text_input("Gmail App Password", type="password")

    st.divider()
    
    # 3. 搜尋設定
    st.subheader("🔍 搜尋條件")
    
    KEYWORDS = {
        "🥚 婦癌 (Gyn Onc)": ["cervical cancer", "ovarian cancer", "endometrial cancer", "immunotherapy", "robotic surgery"],
        "🌊 海扶刀 (HIFU)": ["HIFU", "high intensity focused ultrasound", "uterine leiomyoma", "adenomyosis"],
        "🧬 精準/其他": ["genetic test", "targeted therapy"]
    }
    
    selected_cats = st.multiselect("📚 主題類別", list(KEYWORDS.keys()), default=["🥚 婦癌 (Gyn Onc)"])
    base_keywords = []
    for cat in selected_cats:
        base_keywords.extend(KEYWORDS[cat])
    
    custom_keywords_str = st.text_input("➕ 自訂關鍵字", help="例如: TP53, toxicity")
    if custom_keywords_str:
        base_keywords.extend([k.strip() for k in custom_keywords_str.split(",") if k.strip()])
    
    final_keywords = st.multiselect("最終關鍵字", base_keywords, default=base_keywords)

    use_journals = st.checkbox("限定權威期刊?", value=True)
    PRESET_JOURNALS = ["New England Journal of Medicine", "The Lancet Oncology", "Journal of Clinical Oncology", "Gynecologic Oncology"]
    final_journals = st.multiselect("選擇期刊", PRESET_JOURNALS, default=PRESET_JOURNALS) if use_journals else []

    st.divider()

    # 4. 時間與數量
    date_mode = st.radio("📅 時間模式", ["最近幾天", "指定區間"], index=0)
    date_range_query = ""
    date_params = {}
    
    if date_mode == "最近幾天":
        days_back = st.slider("過去幾天?", 1, 90, 14)
        date_params = {"reldate": days_back}
    else:
        col1, col2 = st.columns(2)
        with col1: day_start = st.number_input("幾天前開始?", 1, 365, 60)
        with col2: day_end = st.number_input("幾天前結束?", 0, 365, 30)
        today = datetime.now()
        d_min = (today - timedelta(days=day_start)).strftime("%Y/%m/%d")
        d_max = (today - timedelta(days=day_end)).strftime("%Y/%m/%d")
        date_range_query = f' AND ("{d_min}"[Date - Publication] : "{d_max}"[Date - Publication])'

    # 數量設定 (0-100)
    max_results = st.number_input("列出篇數上限 (Max Results)", min_value=1, max_value=100, value=20)
    
    # 搜尋按鈕
    if st.button("🚀 極速搜尋 (列出標題)", type="primary", disabled=(not selected_model_name)):
        # 清空舊資料
        st.session_state.articles_data = []
        st.session_state.analysis_cache = {}
        st.session_state.email_queue = []
        st.session_state.search_trigger = True # 觸發搜尋標記

# --- 核心功能函數 ---

def build_query(keywords, journals, date_str_query):
    if not keywords: return ""
    term_q = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    final = term_q
    if journals:
        journal_q = "(" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
        final = f"{term_q} AND {journal_q}"
    if date_str_query: final += date_str_query
    return final

def fetch_headers(query, date_params, limit, email):
    Entrez.email = email
    try:
        search_args = {"db": "pubmed", "term": query, "retmax": limit, "sort": "date"}
        if "reldate" in date_params: search_args["reldate"] = date_params["reldate"]
        
        h = Entrez.esearch(**search_args)
        r = Entrez.read(h)
        ids = r["IdList"]
        if not ids: return []
        
        # 抓取摘要資訊
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
                parsed.append({"id": ids[0], "title":ti, "journal":jo, "abstract":ab, "link":link, "title_zh": "翻譯中..."})
            except: continue
        return parsed
    except Exception as e:
        st.error(f"PubMed Error: {e}"); return []

def batch_translate_titles(articles, key, model_name):
    """
    一次將所有標題打包送給 AI 翻譯 (批次處理，速度極快)
    """
    if not articles: return []
    
    # 準備 Prompt，將標題列表化
    titles_text = "\n".join([f"{i+1}. {art['title']}" for i, art in enumerate(articles)])
    
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    headers = {'Content-Type': 'application/json'}
    
    prompt = f"""
    任務：將以下 {len(articles)} 個醫學論文標題翻譯成「台灣繁體中文」。
    格式：請嚴格按照順序，一行一個翻譯結果，不要有編號，不要有額外文字。
    
    原文標題：
    {titles_text}
    """
    
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            res_text = response.json()['candidates'][0]['content']['parts'][0]['text']
            # 分割回傳的行
            zh_titles = [line.strip() for line in res_text.strip().split('\n') if line.strip()]
            
            # 確保長度一致 (AI有時候會少翻或多編號)
            for i, art in enumerate(articles):
                if i < len(zh_titles):
                    # 移除可能存在的編號 (如 "1. ")
                    clean_title = zh_titles[i].split(". ", 1)[-1] if ". " in zh_titles[i][:4] else zh_titles[i]
                    art['title_zh'] = clean_title
                else:
                    art['title_zh'] = "(翻譯失敗)"
    except:
        pass # 失敗就維持原樣
    return articles

def run_deep_analysis(art, key, model_name):
    """單篇深度分析"""
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    headers = {'Content-Type': 'application/json'}
    
    prompt_text = f"""
    角色：資深婦癌權威醫師。
    任務：針對以下這篇論文進行詳細的學術點評。
    
    標題：{art['title']}
    摘要：{art['abstract']}
    
    請以 HTML 格式 (不含 markdown) 輸出以下結構分析：
    <div style="background:#fff; padding:15px; border:1px solid #ddd; border-radius:5px;">
        <h4 style="color:#2e86c1; margin-top:0;">1. 🧪 研究方法 (Methods)</h4>
        <p>請簡述 Study Design, Patient Population, Intervention。</p>
        
        <h4 style="color:#2e86c1;">2. 💡 發想緣起 (Rationale)</h4>
        <p>推測作者為何進行此研究？解決了什麼臨床痛點？</p>
        
        <h4 style="color:#2e86c1;">3. 📊 結果數據 (Results)</h4>
        <p>請列出關鍵 P-value, HR, OR, Response Rate 等具體數據。</p>
        
        <h4 style="color:#d35400;">4. 🏥 臨床運用與結論 (Clinical Implication)</h4>
        <p>這對婦癌臨床實踐有何具體改變或建議？</p>
    </div>
    """
    
    payload = {"contents": [{"parts": [{"text": prompt_text}]}]}
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            txt = response.json()['candidates'][0]['content']['parts'][0]['text']
            return txt.replace("```html", "").replace("```", "")
        else: return f"<div style='color:red'>分析失敗: {response.text}</div>"
    except Exception as e: return f"<div style='color:red'>錯誤: {str(e)}</div>"

def send_bulk_email(to, pwd, queue):
    if not queue: return False, "清單為空"
    msg = MIMEMultipart()
    msg['From'] = to
    msg['To'] = to
    msg['Subject'] = f"GynOnc 精選文獻彙報 ({len(queue)}篇) - {datetime.now().strftime('%Y-%m-%d')}"
    
    body = "<h2>🧬 您的精選文獻分析報告</h2><hr>"
    for item in queue:
        body += item['html']
        body += "<hr>"
    
    msg.attach(MIMEText(body, 'html'))
    try:
        s = smtplib.SMTP('smtp.gmail.com', 587)
        s.starttls()
        s.login(to, pwd)
        s.send_message(msg); s.quit()
        return True, "已寄出"
    except Exception as e: return False, str(e)

# --- 主程式邏輯 ---

st.title("⚡ GynOnc 文獻快篩系統")
st.caption("先列清單，再點選深入分析")

# 1. 執行搜尋 (只抓標題和摘要，不做分析)
if getattr(st.session_state, 'search_trigger', False):
    with st.status("🔍 正在搜尋 PubMed 並批次翻譯標題...", expanded=True) as status:
        q = build_query(final_keywords, final_journals, date_range_query)
        st.write(f"搜尋語法: `{q[:50]}...`")
        
        # 抓取資料
        raw_articles = fetch_headers(q, date_params, max_results, user_email)
        
        if raw_articles:
            st.write(f"✅ 找到 {len(raw_articles)} 篇，正在進行 AI 標題批次翻譯...")
            # 批次翻譯標題 (這一步很快)
            translated_articles = batch_translate_titles(raw_articles, api_key, selected_model_name)
            st.session_state.articles_data = translated_articles
            status.update(label="搜尋完成！請在下方列表點選查看。", state="complete")
        else:
            status.update(label="❌ 找不到文章", state="error")
    
    st.session_state.search_trigger = False # 關閉觸發器

# 2. 顯示列表 (快篩介面)
if st.session_state.articles_data:
    st.divider()
    st.markdown(f"### 📚 搜尋結果列表 ({len(st.session_state.articles_data)} 篇)")
    
    for i, art in enumerate(st.session_state.articles_data):
        # 使用容器框住每一篇
        with st.container():
            col1, col2 = st.columns([5, 1])
            
            with col1:
                # 顯示標題與中文標題
                st.markdown(f"**{i+1}. {art['title']}**")
                st.markdown(f"<span style='color:#2e86c1; font-size:1.1em;'>{art['title_zh']}</span>", unsafe_allow_html=True)
                st.caption(f"📖 {art['journal']} | 🗓️ [原文連結]({art['link']})")
            
            with col2:
                # 分析按鈕 (Unique key very important)
                btn_key = f"analyze_btn_{i}"
                if st.button("🔍 詳細分析", key=btn_key):
                    # 點擊時，馬上執行分析並存入 Cache
                    with st.spinner("AI 正在深度閱讀此篇文章..."):
                        if art['id'] not in st.session_state.analysis_cache:
                            report = run_deep_analysis(art, api_key, selected_model_name)
                            st.session_state.analysis_cache[art['id']] = report
                            
                            # 自動加入 Email 佇列
                            email_item = {
                                "title": art['title'],
                                "html": f"<h3><a href='{art['link']}'>{art['title']} ({art['title_zh']})</a></h3><p>{art['journal']}</p>{report}"
                            }
                            # 避免重複加入
                            if not any(d['title'] == art['title'] for d in st.session_state.email_queue):
                                st.session_state.email_queue.append(email_item)

            # 如果 Cache 裡有這篇的報告，就展開顯示
            if art['id'] in st.session_state.analysis_cache:
                with st.expander("🩺 查看 AI 深度分析報告", expanded=True):
                    st.markdown(st.session_state.analysis_cache[art['id']], unsafe_allow_html=True)
            
            st.markdown("---")

# 3. 購物車/寄信區
if st.session_state.email_queue:
    st.sidebar.divider()
    st.sidebar.header(f"🛒 已選文獻 ({len(st.session_state.email_queue)})")
    st.sidebar.info("您點擊過「詳細分析」的文章都會自動加入此清單。")
    
    if st.sidebar.button("📩 打包寄出所有已分析文獻", type="primary"):
        if not email_password:
            st.sidebar.error("請輸入 Gmail App Password")
        else:
            ok, msg = send_bulk_email(user_email, email_password, st.session_state.email_queue)
            if ok:
                st.sidebar.success("✅ 郵件已寄出！")
                st.session_state.email_queue = [] # 清空
            else:
                st.sidebar.error(f"❌ 失敗: {msg}")
