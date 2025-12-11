import streamlit as st
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, timedelta
import time
import requests
import json
import concurrent.futures

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 極速文獻系統 v9.0", page_icon="⚡", layout="wide")

# --- Session State ---
if 'articles_data' not in st.session_state: st.session_state.articles_data = []
if 'analysis_cache' not in st.session_state: st.session_state.analysis_cache = {}
if 'email_queue' not in st.session_state: st.session_state.email_queue = []
if 'search_trigger' not in st.session_state: st.session_state.search_trigger = False

# --- 工具函數 ---
def clean_input(text):
    return text.strip() if text else ""

def clean_json_text(text):
    text = text.strip()
    if text.startswith("```json"): text = text[7:]
    elif text.startswith("```"): text = text[3:]
    if text.endswith("```"): text = text[:-3]
    return text.strip()

# --- 側邊欄 ---
with st.sidebar:
    st.header("⚡ 設定與購物車")
    
    # 1. 購物車 (置頂)
    if st.session_state.email_queue:
        with st.expander(f"🛒 購物車 ({len(st.session_state.email_queue)})", expanded=True):
            for item in st.session_state.email_queue:
                st.text(f"• {item['title'][:20]}...")
                
            if 'EMAIL_ADDRESS' in st.secrets: user_email = st.secrets['EMAIL_ADDRESS']
            else: user_email = st.text_input("Email", "lionsmanic@gmail.com")
            
            if 'EMAIL_PASSWORD' in st.secrets: email_password = st.secrets['EMAIL_PASSWORD']
            else: email_password = st.text_input("Gmail App Password", type="password")

            if st.button("📩 寄出", type="primary"):
                st.session_state.trigger_email = True
    else:
        st.info("購物車是空的")
    
    st.divider()

    # 2. API Key
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 Key Ready")
    else:
        raw_key = st.text_input("Gemini API Key", type="password")
        api_key = clean_input(raw_key)

    st.divider()
    
    # 3. 搜尋條件
    st.subheader("🔍 搜尋條件")
    KEYWORDS = {
        "🥚 婦癌": ["cervical cancer", "ovarian cancer", "endometrial cancer", "immunotherapy", "robotic surgery"],
        "🌊 海扶刀": ["HIFU", "high intensity focused ultrasound", "uterine leiomyoma", "adenomyosis"],
        "🧬 精準": ["genetic test", "targeted therapy"]
    }
    
    sel_cat = st.multiselect("類別", list(KEYWORDS.keys()), ["🥚 婦癌"])
    base_k = []
    for c in sel_cat: base_k.extend(KEYWORDS[c])
    cust_k = st.text_input("自訂關鍵字", help="e.g. TP53")
    if cust_k: base_k.extend([k.strip() for k in cust_k.split(",")])
    final_k = st.multiselect("關鍵字", base_k, base_k)

    use_j = st.checkbox("限定期刊", True)
    PRESET_J = ["New England Journal of Medicine", "The Lancet Oncology", "Journal of Clinical Oncology", "Gynecologic Oncology", "Journal of Gynecologic Oncology"]
    final_j = st.multiselect("期刊", PRESET_J, PRESET_J) if use_j else []

    st.divider()

    # 4. 時間設定 (修復部分)
    st.subheader("📅 時間設定")
    date_mode = st.radio("模式", ["最近幾天", "指定區間"], index=0)
    
    date_range_query = ""
    date_params = {} # 用於存放 reldate 等參數

    if date_mode == "最近幾天":
        days_back = st.slider("幾天內?", 1, 90, 14)
        date_params = {"reldate": days_back}
    else:
        col1, col2 = st.columns(2)
        with col1: day_start = st.number_input("幾天前開始?", 1, 365, 60)
        with col2: day_end = st.number_input("幾天前結束?", 0, 365, 30)
        
        # 建立 PubMed 日期區間語法
        today = datetime.now()
        d_min = (today - timedelta(days=day_start)).strftime("%Y/%m/%d")
        d_max = (today - timedelta(days=day_end)).strftime("%Y/%m/%d")
        # 注意：PubMed 語法需要前後空格
        date_range_query = f' AND ("{d_min}"[Date - Publication] : "{d_max}"[Date - Publication])'

    max_res = st.number_input("篇數上限", 1, 100, 20)
    
    if st.button("🚀 極速搜尋", type="primary"):
        if not api_key: st.error("請輸入 API Key")
        else:
            st.session_state.articles_data = []
            st.session_state.search_trigger = True

# --- 核心函數 ---

def build_query(keywords, journals, date_str_query):
    # 基礎關鍵字
    if not keywords: return ""
    term_q = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    
    # 加上期刊
    final = term_q
    if journals:
        journal_q = "(" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
        final = f"{term_q} AND {journal_q}"
    
    # 加上日期區間語法 (如果是指定區間模式)
    if date_str_query:
        final += date_str_query
        
    return final

def fetch_headers(query, date_params, limit, email):
    Entrez.email = email
    try:
        search_args = {"db": "pubmed", "term": query, "retmax": limit, "sort": "date"}
        # 如果是「最近幾天」，加入 reldate 參數
        if "reldate" in date_params: 
            search_args["reldate"] = date_params["reldate"]
        
        # 第一步：搜尋 ID
        h = Entrez.esearch(**search_args)
        r = Entrez.read(h)
        ids = r["IdList"]
        if not ids: return []
        
        # 第二步：抓取詳細資料
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
                link = f"[https://doi.org/](https://doi.org/){doi}" if doi else f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){ids[0]}/"
                parsed.append({"id": ids[0], "title":ti, "journal":jo, "abstract":ab, "link":link, "title_zh": "翻譯中..."})
            except: continue
        return parsed
    except Exception as e:
        st.error(f"PubMed Error: {e}")
        return []

# --- 平行翻譯邏輯 (Multithreading) ---

def translate_chunk(chunk, key):
    if not chunk: return []
    titles = "\n".join([f"{i+1}. {art['title']}" for i, art in enumerate(chunk)])
    url = f"[https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=](https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=){key}"
    headers = {'Content-Type': 'application/json'}
    prompt = f"Translate titles to Traditional Chinese (Taiwan). One per line. No numbering.\n{titles}"
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    try:
        res = requests.post(url, headers=headers, data=json.dumps(payload))
        if res.status_code == 200:
            lines = [l.strip() for l in res.json()['candidates'][0]['content']['parts'][0]['text'].strip().split('\n') if l.strip()]
            for i, art in enumerate(chunk):
                if i < len(lines):
                    clean = lines[i].split(". ", 1)[-1] if ". " in lines[i][:4] else lines[i]
                    art['title_zh'] = clean
                else: art['title_zh'] = art['title']
    except: pass
    return chunk

def batch_translate_parallel(articles, key):
    chunk_size = 10
    chunks = [articles[i:i+chunk_size] for i in range(0, len(articles), chunk_size)]
    results = []
    # 使用 3 個執行緒並發
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(translate_chunk, c, key): c for c in chunks}
        for future in concurrent.futures.as_completed(futures):
            try: results.extend(future.result())
            except: pass
    # 簡單排序回原本順序 (依賴標題匹配)
    title_map = {r['title']: r for r in results}
    final_ordered = []
    for art in articles:
        if art['title'] in title_map: final_ordered.append(title_map[art['title']])
        else: final_ordered.append(art)
    return final_ordered

# --- 深度分析邏輯 (JSON) ---

def run_deep_analysis_json(art, key):
    url = f"[https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=](https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=){key}"
    headers = {'Content-Type': 'application/json'}
    prompt = f"""
    Analyze abstract. Return JSON only. Keys: methods, rationale, results, implication. Values in Traditional Chinese.
    Title: {art['title']}
    Abstract: {art['abstract']}
    """
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    try:
        res = requests.post(url, headers=headers, data=json.dumps(payload))
        if res.status_code == 200:
            txt = clean_json_text(res.json()['candidates'][0]['content']['parts'][0]['text'])
            data = json.loads(txt)
            return f"""
            <div style="background:#fff; padding:15px; border-radius:8px; border:1px solid #eee;">
                <h4 style="color:#2e86c1;">1. 🧪 研究方法</h4><div>{data.get('methods','')}</div>
                <h4 style="color:#2e86c1;">2. 💡 發想緣起</h4><div>{data.get('rationale','')}</div>
                <h4 style="color:#2e86c1;">3. 📊 結果數據</h4><div>{data.get('results','')}</div>
                <h4 style="color:#d35400;">4. 🏥 臨床運用</h4><div>{data.get('implication','')}</div>
            </div>
            """
    except: return "<div style='color:red'>分析失敗</div>"
    return "<div style='color:red'>連線失敗</div>"

def send_mail(to, pwd, queue):
    msg = MIMEMultipart()
    msg['From'] = to
    msg['To'] = to
    msg['Subject'] = f"GynOnc Report {datetime.now().strftime('%Y-%m-%d')}"
    body = "<html><body><h2>文獻報告</h2><hr>" + "".join([i['html'] + "<hr>" for i in queue]) + "</body></html>"
    msg.attach(MIMEText(body, 'html'))
    try:
        s = smtplib.SMTP('smtp.gmail.com', 587); s.starttls()
        s.login(to, pwd); s.send_message(msg); s.quit()
        return True, "OK"
    except Exception as e: return False, str(e)

# --- 主流程 ---

st.title("⚡ GynOnc 極速文獻系統 v9.0")

if st.session_state.search_trigger:
    # 確保 email 變數存在
    search_email = "lionsmanic@gmail.com"
    if 'EMAIL_ADDRESS' in st.secrets: search_email = st.secrets['EMAIL_ADDRESS']
    
    with st.status("🚀 搜尋並平行翻譯中...", expanded=True) as status:
        # 建立查詢
        q = build_query(final_k, final_j, date_range_query)
        st.write(f"語法: `{q[:60]}...`")
        
        # 抓取標題
        raw = fetch_headers(q, date_params, max_res, search_email)
        
        if raw:
            st.write(f"✅ 找到 {len(raw)} 篇，啟動多核心翻譯...")
            # 平行翻譯
            final_list = batch_translate_parallel(raw, api_key)
            st.session_state.articles_data = final_list
            status.update(label="完成！", state="complete")
        else:
            status.update(label="❌ 找不到文章 (請檢查日期或關鍵字)", state="error")
    st.session_state.search_trigger = False

# 顯示結果
if st.session_state.articles_data:
    st.divider()
    for i, art in enumerate(st.session_state.articles_data):
        with st.container():
            c1, c2 = st.columns([5, 1])
            with c1:
                st.markdown(f"**{i+1}. {art['title']}**")
                # 藍色大標題顯示中文
                st.markdown(f"<h4 style='color:#1a5276; margin-top:0;'>{art.get('title_zh', '...')}</h4>", unsafe_allow_html=True)
                st.caption(f"📖 {art['journal']} | [Link]({art['link']})")
            
            with c2:
                if st.button("🔍 詳細分析", key=f"btn_{i}"):
                    with st.spinner("分析中..."):
                        if art['id'] not in st.session_state.analysis_cache:
                            report = run_deep_analysis_json(art, api_key)
                            st.session_state.analysis_cache[art['id']] = report
                            st.session_state.email_queue.append({
                                "title": art['title'],
                                "html": f"<h3>{art['title']}</h3><h4>{art['title_zh']}</h4>{report}"
                            })
                            st.rerun()

            if art['id'] in st.session_state.analysis_cache:
                with st.expander("🩺 深度報告", expanded=True):
                    st.markdown(st.session_state.analysis_cache[art['id']], unsafe_allow_html=True)
            st.markdown("---")

if getattr(st.session_state, 'trigger_email', False):
    m_to = "lionsmanic@gmail.com"
    m_pwd = ""
    if 'EMAIL_ADDRESS' in st.secrets: m_to = st.secrets['EMAIL_ADDRESS']
    if 'EMAIL_PASSWORD' in st.secrets: m_pwd = st.secrets['EMAIL_PASSWORD']
    
    ok, msg = send_mail(m_to, m_pwd, st.session_state.email_queue)
    if ok: 
        st.sidebar.success("已寄出")
        st.session_state.email_queue = []
    else: st.sidebar.error(f"失敗: {msg}")
    st.session_state.trigger_email = False
