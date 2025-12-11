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
st.set_page_config(page_title="GynOnc 文獻系統 v5.0 (穩定版)", page_icon="🛡️", layout="wide")

# --- Session State ---
if 'articles_data' not in st.session_state: st.session_state.articles_data = []
if 'analysis_cache' not in st.session_state: st.session_state.analysis_cache = {}
if 'email_queue' not in st.session_state: st.session_state.email_queue = []
if 'search_trigger' not in st.session_state: st.session_state.search_trigger = False

# --- 核心工具：帶有自動重試功能的 API 呼叫 ---
def call_gemini_api(url, payload, retries=3):
    """
    發送 API 請求，如果遇到 503 (Overloaded) 或 429 (Rate Limit)，
    會自動等待並重試，最多 retries 次。
    """
    headers = {'Content-Type': 'application/json'}
    
    for attempt in range(retries):
        try:
            response = requests.post(url, headers=headers, data=json.dumps(payload))
            
            # 200 OK
            if response.status_code == 200:
                return response.json()
            
            # 503 Service Unavailable (Overloaded) 或 429 Too Many Requests
            elif response.status_code in [503, 429]:
                wait_time = (attempt + 1) * 2  # 第一次等2秒, 第二次等4秒...
                time.sleep(wait_time)
                continue # 重試
            
            # 其他錯誤 (400, 403, 404) -> 直接回傳錯誤，不重試
            else:
                return {"error": f"HTTP {response.status_code}: {response.text}"}
                
        except Exception as e:
            time.sleep(1)
            continue

    return {"error": "Maximum retries exceeded (系統忙碌，請稍後再試)"}

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
    st.header("🛡️ 設定 (穩定版)")
    
    # 1. API Key
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    selected_model_name = None
    if api_key:
        available_models = get_available_models(api_key)
        if available_models:
            default_ix = 0
            if 'gemini-1.5-flash' in available_models: default_ix = available_models.index('gemini-1.5-flash')
            elif 'gemini-pro' in available_models: default_ix = available_models.index('gemini-pro')
            selected_model_name = st.selectbox("✅ AI 模型:", available_models, index=default_ix)

    # 2. Email
    if 'EMAIL_ADDRESS' in st.secrets: user_email = st.secrets['EMAIL_ADDRESS']
    else: user_email = st.text_input("Email", "lionsmanic@gmail.com")
    
    if 'EMAIL_PASSWORD' in st.secrets: email_password = st.secrets['EMAIL_PASSWORD']
    else: email_password = st.text_input("Gmail App Password", type="password")

    st.divider()
    
    # 3. 搜尋條件
    st.subheader("🔍 搜尋條件")
    KEYWORDS = {
        "🥚 婦癌 (Gyn Onc)": ["cervical cancer", "ovarian cancer", "endometrial cancer", "immunotherapy", "robotic surgery"],
        "🌊 海扶刀 (HIFU)": ["HIFU", "high intensity focused ultrasound", "uterine leiomyoma", "adenomyosis"],
        "🧬 精準/其他": ["genetic test", "targeted therapy"]
    }
    
    selected_cats = st.multiselect("📚 主題", list(KEYWORDS.keys()), default=["🥚 婦癌 (Gyn Onc)"])
    base_keywords = []
    for cat in selected_cats: base_keywords.extend(KEYWORDS[cat])
    
    custom_keywords_str = st.text_input("➕ 自訂關鍵字", help="例如: TP53")
    if custom_keywords_str: base_keywords.extend([k.strip() for k in custom_keywords_str.split(",")])
    
    final_keywords = st.multiselect("關鍵字", base_keywords, default=base_keywords)

    use_journals = st.checkbox("限定期刊?", value=True)
    PRESET_JOURNALS = ["New England Journal of Medicine", "The Lancet Oncology", "Journal of Clinical Oncology", "Gynecologic Oncology", "Journal of Gynecologic Oncology"]
    final_journals = st.multiselect("期刊列表", PRESET_JOURNALS, default=PRESET_JOURNALS) if use_journals else []

    st.divider()

    # 4. 時間與數量
    date_mode = st.radio("📅 時間", ["最近幾天", "指定區間"], index=0)
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

    max_results = st.number_input("篇數上限", 1, 100, 20)
    
    if st.button("🚀 極速搜尋", type="primary", disabled=(not selected_model_name)):
        st.session_state.articles_data = []
        st.session_state.analysis_cache = {}
        st.session_state.email_queue = []
        st.session_state.search_trigger = True

# --- 核心函數 ---

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

def chunk_list(lst, n):
    """將列表切分成小塊，避免一次送太多"""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def batch_translate_titles_robust(articles, key, model_name):
    """
    分批翻譯標題，每次只翻譯 5 篇，降低 503 錯誤機率
    """
    if not articles: return []
    
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    
    # 將文章分成每 5 篇一組
    chunk_size = 5
    article_chunks = list(chunk_list(articles, chunk_size))
    
    progress_bar = st.progress(0)
    
    for idx, chunk in enumerate(article_chunks):
        titles_text = "\n".join([f"{i+1}. {art['title']}" for i, art in enumerate(chunk)])
        
        prompt = f"""
        任務：翻譯以下 {len(chunk)} 個醫學標題為繁體中文。
        格式：一行一個結果，嚴禁編號，嚴禁多餘文字。
        原文：
        {titles_text}
        """
        payload = {"contents": [{"parts": [{"text": prompt}]}]}
        
        # 使用帶有重試機制的呼叫
        result = call_gemini_api(url, payload)
        
        if "error" not in result:
            try:
                res_text = result['candidates'][0]['content']['parts'][0]['text']
                zh_titles = [line.strip() for line in res_text.strip().split('\n') if line.strip()]
                
                for i, art in enumerate(chunk):
                    if i < len(zh_titles):
                        clean = zh_titles[i].split(". ", 1)[-1] if ". " in zh_titles[i][:4] else zh_titles[i]
                        art['title_zh'] = clean
                    else:
                        art['title_zh'] = "(翻譯格式錯誤)"
            except:
                for art in chunk: art['title_zh'] = "(解析失敗)"
        else:
            for art in chunk: art['title_zh'] = "(翻譯連線逾時)"
            
        # 更新進度條
        progress_bar.progress((idx + 1) / len(article_chunks))
        time.sleep(0.5) # 稍微休息一下，對 API 溫柔一點
        
    return articles

def run_deep_analysis_robust(art, key, model_name):
    """
    深度分析 (含重試機制 + HTML 強制格式)
    """
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    
    prompt_text = f"""
    角色：資深婦癌權威醫師。
    標題：{art['title']}
    摘要：{art['abstract']}
    
    【輸出格式要求】：
    1. 輸出 **純 HTML**。
    2. **嚴禁** Markdown。
    3. 所有標題用 <h4 style="color:#2e86c1;">。
    4. 全部包在 <div> 內。
    
    模板：
    <div style="font-family: sans-serif; line-height: 1.6;">
        <h4 style="color:#2e86c1; margin-top:0; border-bottom: 2px solid #eee;">1. 🧪 研究方法 (Methods)</h4>
        <p>簡述 Study Design, Patient Population。</p>
        
        <h4 style="color:#2e86c1; border-bottom: 2px solid #eee;">2. 💡 發想緣起 (Rationale)</h4>
        <p>為何做此研究？解決什麼痛點？</p>
        
        <h4 style="color:#2e86c1; border-bottom: 2px solid #eee;">3. 📊 結果數據 (Results)</h4>
        <ul><li>關鍵數據 (P-value, HR)...</li></ul>
        
        <h4 style="color:#d35400; border-bottom: 2px solid #eee;">4. 🏥 臨床運用 (Implication)</h4>
        <p>臨床建議。</p>
    </div>
    """
    
    payload = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    # 呼叫 API (帶重試)
    result = call_gemini_api(url, payload)
    
    if "error" in result:
        # 如果重試多次還是失敗，回傳紅字錯誤
        return f"<div style='color:red; border:1px solid red; padding:10px;'>❌ 分析失敗 (系統忙碌): {result['error']}</div>"
    
    try:
        txt = result['candidates'][0]['content']['parts'][0]['text']
        return txt.replace("```html", "").replace("```", "").strip()
    except Exception as e:
        return f"<div style='color:red'>解析錯誤: {str(e)}</div>"

def send_bulk_email(to, pwd, queue):
    if not queue: return False, "清單為空"
    msg = MIMEMultipart()
    msg['From'] = to
    msg['To'] = to
    msg['Subject'] = f"GynOnc 文獻彙報 ({len(queue)}篇) - {datetime.now().strftime('%Y-%m-%d')}"
    
    body = """
    <html><body style="font-family: Arial, sans-serif; color: #333;">
    <h2 style="color: #2c3e50;">🧬 文獻分析報告</h2>
    <hr>
    """
    for item in queue:
        body += item['html']
        body += "<hr style='margin: 30px 0; border: 0; border-top: 1px solid #eee;'>"
    body += "</body></html>"
    
    msg.attach(MIMEText(body, 'html'))
    try:
        s = smtplib.SMTP('smtp.gmail.com', 587)
        s.starttls()
        s.login(to, pwd)
        s.send_message(msg); s.quit()
        return True, "已寄出"
    except Exception as e: return False, str(e)

# --- 主程式 ---

st.title("🛡️ GynOnc 文獻系統 v5.0")
st.caption("穩定版：內建 503 自動重試與分批處理機制")

# 1. 搜尋
if st.session_state.search_trigger:
    with st.status("🔍 正在執行穩定搜尋...", expanded=True) as status:
        q = build_query(final_keywords, final_journals, date_range_query)
        st.write(f"搜尋語法: `{q[:50]}...`")
        
        raw_articles = fetch_headers(q, date_params, max_results, user_email)
        
        if raw_articles:
            st.write(f"✅ 找到 {len(raw_articles)} 篇")
            st.write("🔄 正在分批翻譯標題 (每5篇一組，避免當機)...")
            
            # 使用分批翻譯函數
            translated_articles = batch_translate_titles_robust(raw_articles, api_key, selected_model_name)
            
            st.session_state.articles_data = translated_articles
            status.update(label="搜尋完成！", state="complete")
        else:
            status.update(label="❌ 找不到文章", state="error")
    
    st.session_state.search_trigger = False

# 2. 列表
if st.session_state.articles_data:
    st.divider()
    st.markdown(f"### 📚 結果列表 ({len(st.session_state.articles_data)} 篇)")
    
    for i, art in enumerate(st.session_state.articles_data):
        with st.container():
            col1, col2 = st.columns([5, 1])
            with col1:
                st.markdown(f"**{i+1}. {art['title']}**")
                # 顯示中文標題，若失敗會顯示原因
                st.markdown(f"<span style='color:#2e86c1; font-size:1.1em;'>{art.get('title_zh', '翻譯中...')}</span>", unsafe_allow_html=True)
                st.caption(f"📖 {art['journal']} | [連結]({art['link']})")
            
            with col2:
                btn_key = f"analyze_{art['id']}_{i}"
                if st.button("🔍 詳細分析", key=btn_key):
                    # 執行分析 (帶有重試機制)
                    with st.spinner("AI 正在深度閱讀 (若忙碌將自動重試)..."):
                        if art['id'] not in st.session_state.analysis_cache:
                            report = run_deep_analysis_robust(art, api_key, selected_model_name)
                            st.session_state.analysis_cache[art['id']] = report
                            
                            email_item = {
                                "title": art['title'],
                                "html": f"""
                                <div style="background-color: #f9f9f9; padding: 20px; border-radius: 5px; margin-bottom: 20px;">
                                    <h3 style="margin-top: 0; color: #1a5276;"><a href='{art['link']}' style="text-decoration: none; color: #1a5276;">{art['title']}</a></h3>
                                    <h4 style="margin-top: 5px; color: #2e86c1;">{art.get('title_zh', '')}</h4>
                                    <p style="color: #666; font-size: 0.9em;">📖 {art['journal']}</p>
                                    {report}
                                </div>
                                """
                            }
                            if not any(d['title'] == art['title'] for d in st.session_state.email_queue):
                                st.session_state.email_queue.append(email_item)

            if art['id'] in st.session_state.analysis_cache:
                with st.expander("🩺 分析報告", expanded=True):
                    st.markdown(st.session_state.analysis_cache[art['id']], unsafe_allow_html=True)
            st.markdown("---")

# 3. 寄信
if st.session_state.email_queue:
    st.sidebar.divider()
    st.sidebar.header(f"🛒 購物車 ({len(st.session_state.email_queue)})")
    if st.sidebar.button("📩 打包寄出"):
        if not email_password: st.sidebar.error("缺密碼")
        else:
            ok, msg = send_bulk_email(user_email, email_password, st.session_state.email_queue)
            if ok: 
                st.sidebar.success("已寄出！")
                st.session_state.email_queue = []
            else: st.sidebar.error(msg)
