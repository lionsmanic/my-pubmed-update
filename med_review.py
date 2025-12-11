import streamlit as st
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, timedelta
import time
import requests
import json
import concurrent.futures # 引入平行處理庫

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 極速翻譯版 v8.0", page_icon="⚡", layout="wide")

# --- Session State ---
if 'articles_data' not in st.session_state: st.session_state.articles_data = []
if 'analysis_cache' not in st.session_state: st.session_state.analysis_cache = {}
if 'email_queue' not in st.session_state: st.session_state.email_queue = []
if 'search_trigger' not in st.session_state: st.session_state.search_trigger = False

# --- 工具函數 ---

def clean_input(text):
    """清理 API Key，去除空格"""
    return text.strip() if text else ""

def clean_json_text(text):
    """清理 JSON 格式"""
    text = text.strip()
    if text.startswith("```json"): text = text[7:]
    elif text.startswith("```"): text = text[3:]
    if text.endswith("```"): text = text[:-3]
    return text.strip()

# --- 核心邏輯：多執行緒翻譯 ---

def translate_chunk(chunk, key):
    """
    翻譯一個小區塊 (Worker Function)
    """
    if not chunk: return []
    
    # 組合 Prompt
    titles_text = "\n".join([f"{i+1}. {art['title']}" for i, art in enumerate(chunk)])
    url = f"[https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=](https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=){key}"
    headers = {'Content-Type': 'application/json'}
    
    prompt = f"""
    Translate these medical titles to Traditional Chinese (Taiwan).
    Format: One translation per line. No numbering. No English.
    
    Titles:
    {titles_text}
    """
    
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            res_text = response.json()['candidates'][0]['content']['parts'][0]['text']
            # 分割回傳的每一行
            lines = [l.strip() for l in res_text.strip().split('\n') if l.strip()]
            
            # 對應回文章
            for i, art in enumerate(chunk):
                if i < len(lines):
                    # 移除可能被 AI 加上的編號
                    clean_zh = lines[i].split(". ", 1)[-1] if ". " in lines[i][:4] else lines[i]
                    art['title_zh'] = clean_zh
                else:
                    art['title_zh'] = art['title'] # 沒翻到就顯示原文
        else:
            for art in chunk: art['title_zh'] = "(翻譯忙碌中)"
    except:
        for art in chunk: art['title_zh'] = "(連線錯誤)"
        
    return chunk

def batch_translate_parallel(articles, key):
    """
    主控台：將文章分組，並發送給多個 Worker 同時跑
    """
    chunk_size = 10 # 每一組 10 篇
    chunks = [articles[i:i + chunk_size] for i in range(0, len(articles), chunk_size)]
    
    results = []
    
    # 開啟 3 個執行緒 (Thread) 同時跑
    # 注意：免費版 API 限制較多，開太多會被擋，3 個是安全值
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        # 提交任務
        future_to_chunk = {executor.submit(translate_chunk, chunk, key): chunk for chunk in chunks}
        
        # 等待結果
        for future in concurrent.futures.as_completed(future_to_chunk):
            try:
                data = future.result()
                results.extend(data)
            except Exception as e:
                pass
                
    # 因為多執行緒回傳順序不固定，這裡簡單處理直接回傳結果列表
    # (如果很在意順序，可以用 map，但 list extend 夠用了)
    return sorted(results, key=lambda x: articles.index(x)) # 嘗試恢復原本順序 (依賴 object id 或內容)

# --- 側邊欄 ---
with st.sidebar:
    st.header("⚡ 極速翻譯設定")
    
    # 1. 購物車
    if st.session_state.email_queue:
        with st.expander(f"🛒 購物車 ({len(st.session_state.email_queue)})", expanded=True):
            if 'EMAIL_ADDRESS' in st.secrets: user_email = st.secrets['EMAIL_ADDRESS']
            else: user_email = st.text_input("Email", "lionsmanic@gmail.com")
            
            if 'EMAIL_PASSWORD' in st.secrets: email_password = st.secrets['EMAIL_PASSWORD']
            else: email_password = st.text_input("Gmail App Password", type="password")

            if st.button("📩 寄出", type="primary"):
                st.session_state.trigger_email = True
    
    st.divider()

    # 2. API Key (Auto Clean)
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 Key Ready")
    else:
        raw_key = st.text_input("Gemini API Key", type="password")
        api_key = clean_input(raw_key)

    st.divider()
    
    # 3. 搜尋
    st.subheader("🔍 條件")
    KEYWORDS = {
        "🥚 婦癌": ["cervical cancer", "ovarian cancer", "endometrial cancer", "immunotherapy", "robotic surgery"],
        "🌊 海扶刀": ["HIFU", "high intensity focused ultrasound", "uterine leiomyoma", "adenomyosis"],
        "🧬 精準": ["genetic test", "targeted therapy"]
    }
    
    sel_cat = st.multiselect("類別", list(KEYWORDS.keys()), ["🥚 婦癌"])
    base_k = []
    for c in sel_cat: base_k.extend(KEYWORDS[c])
    cust_k = st.text_input("自訂", help="e.g. TP53")
    if cust_k: base_k.extend([k.strip() for k in cust_k.split(",")])
    final_k = st.multiselect("關鍵字", base_k, base_k)

    use_j = st.checkbox("限定期刊", True)
    PRESET_J = ["New England Journal of Medicine", "The Lancet Oncology", "Journal of Clinical Oncology", "Gynecologic Oncology", "Journal of Gynecologic Oncology"]
    final_j = st.multiselect("期刊", PRESET_J, PRESET_J) if use_j else []

    days_back = st.slider("幾天內?", 1, 90, 14)
    max_res = st.number_input("篇數", 1, 50, 10) # 建議 20 篇內體驗最好
    
    if st.button("🚀 搜尋並翻譯", type="primary"):
        if not api_key: st.error("缺 API Key")
        else:
            st.session_state.articles_data = []
            st.session_state.search_trigger = True

# --- 主程式函數 ---

def build_query(keywords, journals, days):
    q = "(" + " OR ".join([f'"{k}"[Title/Abstract]' for k in keywords]) + ")"
    if journals: q = f"{q} AND (" + " OR ".join([f'"{j}"[Journal]' for j in journals]) + ")"
    return q

def fetch_headers(query, days, limit, email):
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
                link = f"[https://doi.org/](https://doi.org/){doi}" if doi else f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){ids[0]}/"
                parsed.append({"id": ids[0], "title":ti, "journal":jo, "abstract":ab, "link":link, "title_zh": "翻譯中..."})
            except: continue
        return parsed
    except: return []

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

st.title("⚡ GynOnc 極速翻譯版 v8.0")

if st.session_state.search_trigger:
    # 為了避免變數錯誤，這裡定義臨時 email
    temp_email = "lionsmanic@gmail.com"
    if 'EMAIL_ADDRESS' in st.secrets: temp_email = st.secrets['EMAIL_ADDRESS']
    
    with st.status("🚀 啟動多核心引擎：搜尋 + 平行翻譯中...", expanded=True) as status:
        q = build_query(final_k, final_j, days_back)
        raw = fetch_headers(q, {"reldate": days_back}, max_res, temp_email)
        
        if raw:
            st.write(f"✅ 找到 {len(raw)} 篇，正在同時翻譯...")
            # 關鍵：平行翻譯
            final_list = batch_translate_parallel(raw, api_key)
            st.session_state.articles_data = final_list
            status.update(label="完成！", state="complete")
        else:
            status.update(label="無結果", state="error")
    st.session_state.search_trigger = False

# 顯示列表
if st.session_state.articles_data:
    st.divider()
    for i, art in enumerate(st.session_state.articles_data):
        with st.container():
            c1, c2 = st.columns([5, 1])
            with c1:
                # 這裡直接顯示翻譯好的標題
                st.markdown(f"**{i+1}. {art['title']}**")
                # 藍色大字體顯示中文標題
                st.markdown(f"<h4 style='color:#1a5276; margin-top:0;'>{art.get('title_zh', '...')}</h4>", unsafe_allow_html=True)
                st.caption(f"📖 {art['journal']} | [Link]({art['link']})")
            
            with c2:
                if st.button("🔍 詳細分析", key=f"btn_{i}"):
                    with st.spinner("分析中..."):
                        if art['id'] not in st.session_state.analysis_cache:
                            report = run_deep_analysis_json(art, api_key)
                            st.session_state.analysis_cache[art['id']] = report
                            
                            # 加入購物車
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
    # 再次獲取密碼 (需要確保 sidebar 變數可及性，或使用 secrets)
    m_to = "lionsmanic@gmail.com"
    m_pwd = ""
    if 'EMAIL_ADDRESS' in st.secrets: m_to = st.secrets['EMAIL_ADDRESS']
    if 'EMAIL_PASSWORD' in st.secrets: m_pwd = st.secrets['EMAIL_PASSWORD']
    
    # 如果沒設定 secrets，這裡會失敗，建議正式環境務必設 secrets
    ok, msg = send_mail(m_to, m_pwd, st.session_state.email_queue)
    if ok: 
        st.sidebar.success("已寄出")
        st.session_state.email_queue = []
    else: st.sidebar.error(f"失敗: {msg}")
    st.session_state.trigger_email = False
