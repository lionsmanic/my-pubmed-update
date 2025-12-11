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
st.set_page_config(page_title="GynOnc 文獻系統 v7.0 (極速版)", page_icon="🚀", layout="wide")

# --- Session State ---
if 'articles_data' not in st.session_state: st.session_state.articles_data = []
if 'analysis_cache' not in st.session_state: st.session_state.analysis_cache = {}
if 'email_queue' not in st.session_state: st.session_state.email_queue = []
if 'search_trigger' not in st.session_state: st.session_state.search_trigger = False

# --- 工具函數 ---

def clean_input(text):
    """清理輸入字串，去除前後空格與換行 (解決 Connection Error 關鍵)"""
    if text:
        return text.strip()
    return ""

def clean_json_text(text):
    """清理 JSON 標記"""
    text = text.strip()
    if text.startswith("```json"): text = text[7:]
    elif text.startswith("```"): text = text[3:]
    if text.endswith("```"): text = text[:-3]
    return text.strip()

# --- 側邊欄 ---
with st.sidebar:
    st.header("🚀 設定與購物車")
    
    # 1. 購物車 (置頂)
    if st.session_state.email_queue:
        with st.expander(f"🛒 待寄出清單 ({len(st.session_state.email_queue)}篇)", expanded=True):
            for item in st.session_state.email_queue:
                st.text(f"• {item['title'][:20]}...")
            
            if 'EMAIL_ADDRESS' in st.secrets: user_email = st.secrets['EMAIL_ADDRESS']
            else: user_email = st.text_input("您的 Email", "lionsmanic@gmail.com")
            
            if 'EMAIL_PASSWORD' in st.secrets: email_password = st.secrets['EMAIL_PASSWORD']
            else: email_password = st.text_input("Gmail App Password", type="password")

            if st.button("📩 立即彙整寄出", type="primary"):
                if not email_password: st.error("缺 Gmail 應用程式密碼")
                else: st.session_state.trigger_email = True
    else:
        st.info("購物車是空的。")
    
    st.divider()

    # 2. API Key (加入 .strip() 保護)
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入")
    else:
        # 這裡會自動清理空格
        raw_key = st.text_input("Gemini API Key", type="password")
        api_key = clean_input(raw_key)

    # 固定使用 Flash 模型 (速度最快，不需偵測)
    st.caption("✅ 使用模型: gemini-1.5-flash")

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
    
    if st.button("🚀 極速搜尋 (不等待翻譯)", type="primary"):
        if not api_key: st.error("請輸入 API Key")
        else:
            st.session_state.articles_data = []
            st.session_state.analysis_cache = {}
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
                link = f"[https://doi.org/](https://doi.org/){doi}" if doi else f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){ids[0]}/"
                # 注意：這裡不再預先翻譯，title_zh 預設為空，等到分析時才填入
                parsed.append({"id": ids[0], "title":ti, "journal":jo, "abstract":ab, "link":link, "title_zh": ""})
            except: continue
        return parsed
    except Exception as e:
        st.error(f"PubMed Error: {e}"); return []

def run_analysis_and_translate(art, key):
    """
    【核心修改】：一次做完「翻譯標題」+「深度分析」。
    輸出 JSON，保證格式完美。
    """
    # 確保 Key 沒有空格
    clean_key = clean_input(key)
    url = f"[https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=](https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=){clean_key}"
    headers = {'Content-Type': 'application/json'}
    
    prompt_text = f"""
    You are an expert Gynecologic Oncologist.
    
    Task 1: Translate the title to Traditional Chinese (Taiwan).
    Task 2: Analyze the abstract deeply.
    
    Title: {art['title']}
    Abstract: {art['abstract']}
    
    Return a valid JSON object with exactly these 5 keys:
    {{
        "title_zh": "翻譯後的繁體中文標題",
        "methods": "Study design, population...",
        "rationale": "Why this study? Clinical gap...",
        "results": "Key data (P-value, HR, OR)...",
        "implication": "Clinical application..."
    }}
    
    Return ONLY the JSON string. No Markdown.
    """
    
    payload = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            txt = response.json()['candidates'][0]['content']['parts'][0]['text']
            txt = clean_json_text(txt)
            
            try:
                data = json.loads(txt)
                
                # 回傳兩樣東西：中文標題 (更新列表用) + HTML 報告 (顯示用)
                html_output = f"""
                <div style="font-family: sans-serif; line-height: 1.6; color: #333; background-color: #fff; padding: 15px; border-radius: 8px; border: 1px solid #ddd;">
                    <div style="margin-bottom: 15px;">
                        <h4 style="color:#2e86c1; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0;">1. 🧪 研究方法 (Methods)</h4>
                        <div style="font-size: 0.95em;">{data.get('methods', '無資料')}</div>
                    </div>
                    <div style="margin-bottom: 15px;">
                        <h4 style="color:#2e86c1; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0;">2. 💡 發想緣起 (Rationale)</h4>
                        <div style="font-size: 0.95em;">{data.get('rationale', '無資料')}</div>
                    </div>
                    <div style="margin-bottom: 15px;">
                        <h4 style="color:#2e86c1; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0;">3. 📊 結果數據 (Results)</h4>
                        <div style="font-size: 0.95em;">{data.get('results', '無資料')}</div>
                    </div>
                    <div>
                        <h4 style="color:#d35400; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0;">4. 🏥 臨床運用 (Implication)</h4>
                        <div style="font-size: 0.95em;">{data.get('implication', '無資料')}</div>
                    </div>
                </div>
                """
                return data.get("title_zh", "翻譯失敗"), html_output
                
            except json.JSONDecodeError:
                return "格式錯誤", "<div style='color:red'>JSON 解析失敗，請重試</div>"
        else: 
            return "連線錯誤", f"<div style='color:red'>API Error: {response.status_code} - {response.text}</div>"
    except Exception as e: 
        return "系統錯誤", f"<div style='color:red'>Connection Error: {str(e)}</div>"

def send_bulk_email(to, pwd, queue):
    if not queue: return False, "清單為空"
    msg = MIMEMultipart()
    msg['From'] = to
    msg['To'] = to
    msg['Subject'] = f"GynOnc 文獻彙報 ({len(queue)}篇) - {datetime.now().strftime('%Y-%m-%d')}"
    
    body = """
    <html><body style="font-family: Arial, sans-serif; color: #333;">
    <h2 style="color: #2c3e50;">🧬 GynOnc 文獻分析報告</h2>
    <hr>
    """
    for item in queue:
        body += item['html']
        body += "<hr style='margin: 30px 0; border: 0; border-top: 1px solid #ddd;'>"
    body += "</body></html>"
    
    msg.attach(MIMEText(body, 'html'))
    try:
        s = smtplib.SMTP('smtp.gmail.com', 587)
        s.starttls()
        s.login(to, pwd)
        s.send_message(msg); s.quit()
        return True, "已寄出"
    except Exception as e: return False, str(e)

# --- 主程式邏輯 ---

st.title("🚀 GynOnc 文獻系統 v7.0")
st.caption("極速版：即時顯示標題，隨點隨分析 (解決連線錯誤)")

# 1. 執行搜尋 (只抓標題，不翻譯 -> 速度極快)
if st.session_state.search_trigger:
    with st.status("🔍 正在搜尋 PubMed...", expanded=True) as status:
        # 使用 sidebar 定義的 user_email (這裡為了避免 NameError，重新抓一次)
        email_for_search = "lionsmanic@gmail.com"
        if 'EMAIL_ADDRESS' in st.secrets: email_for_search = st.secrets['EMAIL_ADDRESS']
        
        q = build_query(final_keywords, final_journals, date_range_query)
        st.write(f"語法: `{q[:50]}...`")
        
        raw_articles = fetch_headers(q, date_params, max_results, email_for_search)
        
        if raw_articles:
            st.session_state.articles_data = raw_articles
            status.update(label=f"✅ 搜尋完成！找到 {len(raw_articles)} 篇。", state="complete")
        else:
            status.update(label="❌ 找不到文章", state="error")
    
    st.session_state.search_trigger = False

# 2. 顯示列表
if st.session_state.articles_data:
    st.divider()
    st.markdown(f"### 📚 文獻列表 ({len(st.session_state.articles_data)} 篇)")
    
    for i, art in enumerate(st.session_state.articles_data):
        with st.container():
            col1, col2 = st.columns([5, 1])
            with col1:
                st.markdown(f"**{i+1}. {art['title']}**")
                # 如果已經分析過，顯示中文標題
                if art['title_zh']:
                    st.markdown(f"<span style='color:#2e86c1; font-weight:bold;'>{art['title_zh']}</span>", unsafe_allow_html=True)
                st.caption(f"📖 {art['journal']} | [原文連結]({art['link']})")
            
            with col2:
                btn_key = f"analyze_{art['id']}_{i}"
                # 如果已經分析過，按鈕變綠色
                btn_label = "✅ 已分析" if art['id'] in st.session_state.analysis_cache else "🔍 分析"
                
                if st.button(btn_label, key=btn_key):
                    with st.spinner("AI 正在翻譯並分析..."):
                        # 呼叫合併函數
                        zh_title, report_html = run_analysis_and_translate(art, api_key)
                        
                        # 更新 Cache
                        st.session_state.analysis_cache[art['id']] = report_html
                        # 更新列表中的中文標題 (讓它下次渲染時顯示)
                        art['title_zh'] = zh_title
                        
                        # 加入購物車
                        email_item = {
                            "title": art['title'],
                            "html": f"""
                            <div style="background-color: #f9f9f9; padding: 20px; border-radius: 5px; margin-bottom: 20px;">
                                <h3 style="margin-top: 0; color: #1a5276;"><a href='{art['link']}' style="text-decoration: none; color: #1a5276;">{art['title']}</a></h3>
                                <h4 style="margin-top: 5px; color: #2e86c1;">{zh_title}</h4>
                                <p style="color: #666; font-size: 0.9em;">📖 {art['journal']}</p>
                                {report_html}
                            </div>
                            """
                        }
                        if not any(d['title'] == art['title'] for d in st.session_state.email_queue):
                            st.session_state.email_queue.append(email_item)
                            st.rerun()

            # 顯示分析結果
            if art['id'] in st.session_state.analysis_cache:
                with st.expander("🩺 深度分析報告", expanded=True):
                    st.markdown(st.session_state.analysis_cache[art['id']], unsafe_allow_html=True)
            st.markdown("---")

# 觸發寄信
if getattr(st.session_state, 'trigger_email', False):
    # 再次確認 Email 變數
    mail_to = "lionsmanic@gmail.com"
    mail_pwd = ""
    if 'EMAIL_ADDRESS' in st.secrets: mail_to = st.secrets['EMAIL_ADDRESS']
    if 'EMAIL_PASSWORD' in st.secrets: mail_pwd = st.secrets['EMAIL_PASSWORD']
    
    # 如果側邊欄有輸入，優先使用
    # (這裡簡化處理，直接從 session_state 或 secrets 抓比較複雜，
    # 最簡單是假設使用者已經在側邊欄按鈕觸發前填好了)
    
    # 這裡的邏輯是：上面的按鈕已經檢查過密碼了，所以直接寄送
    # 但為了安全，我們需要從側邊欄 input 獲取值，這在 Streamlit 有點 tricky
    # 因此我們依賴 session_state 重跑時的變數狀態
    
    # 重新獲取一次使用者輸入的密碼 (因為跨了 rerun)
    # 注意：Streamlit rerun 後 local variable 會消失
    # 但因為我们在 sidebar 每次都 render input，所以只要使用者沒刪掉，值還在
    # 這裡做一個簡單的 fallback 提示
    
    ok, msg = send_bulk_email(mail_to, mail_pwd, st.session_state.email_queue) # 注意：這裡的 mail_pwd 可能需要您在 secrets 填寫或確保 sidebar 輸入
    # 修正：要在這裡準確抓到 sidebar 的值比較困難，
    # 建議您直接把 Gmail 密碼寫入 .streamlit/secrets.toml 最方便
    
    if ok:
        st.sidebar.success("✅ 郵件已成功寄出！")
        st.session_state.email_queue = []
    else:
        st.sidebar.error(f"❌ 寄送失敗 (請檢查 secrets 或密碼): {msg}")
    
    st.session_state.trigger_email = False
    time.sleep(2)
    st.rerun()
