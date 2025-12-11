import streamlit as st
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, timedelta
import time
import requests
import json
import re

# --- 頁面設定 ---
st.set_page_config(page_title="GynOnc 文獻系統 v6.2 (修復版)", page_icon="💎", layout="wide")

# --- Session State 初始化 ---
if 'articles_data' not in st.session_state: st.session_state.articles_data = []
if 'analysis_cache' not in st.session_state: st.session_state.analysis_cache = {}
if 'email_queue' not in st.session_state: st.session_state.email_queue = []
if 'search_trigger' not in st.session_state: st.session_state.search_trigger = False

# --- 核心工具函數 ---

def clean_json_text(text):
    """清理 AI 回傳的 JSON 字串"""
    text = text.strip()
    if text.startswith("```json"): text = text[7:]
    elif text.startswith("```"): text = text[3:]
    if text.endswith("```"): text = text[:-3]
    return text.strip()

def get_available_models(api_key):
    """嘗試取得可用模型列表"""
    url = f"[https://generativelanguage.googleapis.com/v1beta/models?key=](https://generativelanguage.googleapis.com/v1beta/models?key=){api_key}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            models = [m['name'].replace('models/', '') for m in data.get('models', []) 
                      if 'generateContent' in m.get('supportedGenerationMethods', [])]
            return models
        return []
    except: return []

# --- 側邊欄設定區 ---
with st.sidebar:
    st.header("💎 設定與購物車")
    
    # 1. Email 設定 (移到最外層，確保變數一定存在)
    if 'EMAIL_ADDRESS' in st.secrets: 
        user_email = st.secrets['EMAIL_ADDRESS']
    else: 
        user_email = st.text_input("您的 Email (必填)", "lionsmanic@gmail.com")
        
    if 'EMAIL_PASSWORD' in st.secrets: 
        email_password = st.secrets['EMAIL_PASSWORD']
    else: 
        email_password = st.text_input("Gmail App Password (寄信用)", type="password")
    
    st.divider()

    # 2. 購物車顯示區
    if st.session_state.email_queue:
        with st.expander(f"🛒 待寄出清單 ({len(st.session_state.email_queue)}篇)", expanded=True):
            for item in st.session_state.email_queue:
                st.text(f"• {item['title'][:20]}...")
            
            # 觸發寄信
            if st.button("📩 立即彙整寄出", type="primary"):
                if not email_password:
                    st.error("請輸入 Gmail 應用程式密碼")
                else:
                    st.session_state.trigger_email = True
    else:
        st.info("目前購物車是空的。請在右側點擊「詳細分析」加入文章。")

    st.divider()

    # 3. API Key 與模型
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    selected_model_name = "gemini-1.5-flash" # 預設值

    if api_key:
        with st.spinner("連線 Google 中..."):
            available_models = get_available_models(api_key)
        
        if available_models:
            st.success("✅ 連線成功")
            default_ix = 0
            if 'gemini-1.5-flash' in available_models: default_ix = available_models.index('gemini-1.5-flash')
            elif 'gemini-pro' in available_models: default_ix = available_models.index('gemini-pro')
            selected_model_name = st.selectbox("選擇模型:", available_models, index=default_ix)
        else:
            st.warning("⚠️ 無法自動取得清單，已使用預設值。")

    st.divider()
    
    # 4. 搜尋條件
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

    # 5. 時間與數量
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
    
    if st.button("🚀 極速搜尋", type="primary"):
        if not api_key:
            st.error("請先輸入 API Key")
        else:
            st.session_state.articles_data = []
            st.session_state.analysis_cache = {}
            st.session_state.search_trigger = True

# --- 主要功能函數 ---

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
                parsed.append({"id": ids[0], "title":ti, "journal":jo, "abstract":ab, "link":link, "title_zh": "翻譯中..."})
            except: continue
        return parsed
    except Exception as e:
        st.error(f"PubMed Error: {e}"); return []

def chunk_list(lst, n):
    for i in range(0, len(lst), n): yield lst[i:i + n]

def batch_translate_titles_robust(articles, key, model_name):
    if not articles: return []
    url = f"[https://generativelanguage.googleapis.com/v1beta/models/](https://generativelanguage.googleapis.com/v1beta/models/){model_name}:generateContent?key={key}"
    headers = {'Content-Type': 'application/json'}
    
    chunk_size = 5
    article_chunks = list(chunk_list(articles, chunk_size))
    progress_bar = st.progress(0)
    
    for idx, chunk in enumerate(article_chunks):
        titles_text = "\n".join([f"{i+1}. {art['title']}" for i, art in enumerate(chunk)])
        prompt = f"任務：翻譯醫學標題為繁體中文。\n格式：一行一個結果，無編號。\n原文：\n{titles_text}"
        payload = {"contents": [{"parts": [{"text": prompt}]}]}
        
        try:
            response = requests.post(url, headers=headers, data=json.dumps(payload))
            if response.status_code == 200:
                res_text = response.json()['candidates'][0]['content']['parts'][0]['text']
                zh_titles = [line.strip() for line in res_text.strip().split('\n') if line.strip()]
                for i, art in enumerate(chunk):
                    if i < len(zh_titles):
                        clean = zh_titles[i].split(". ", 1)[-1] if ". " in zh_titles[i][:4] else zh_titles[i]
                        art['title_zh'] = clean
                    else: art['title_zh'] = art['title']
        except: pass
        
        progress_bar.progress((idx + 1) / len(article_chunks))
        time.sleep(0.5)
        
    return articles

def run_deep_analysis_json(art, key, model_name):
    """AI 輸出 JSON -> Python 轉成 HTML (解決亂碼)"""
    url = f"[https://generativelanguage.googleapis.com/v1beta/models/](https://generativelanguage.googleapis.com/v1beta/models/){model_name}:generateContent?key={key}"
    headers = {'Content-Type': 'application/json'}
    
    prompt_text = f"""
    You are an expert Gynecologic Oncologist. Analyze this abstract.
    Title: {art['title']}
    Abstract: {art['abstract']}
    
    Return a valid JSON object with exactly these 4 keys (value must be Traditional Chinese string):
    {{
        "methods": "簡述研究設計、收案對象...",
        "rationale": "發想緣起、臨床痛點...",
        "results": "列點說明關鍵數據 (P值, HR)...",
        "implication": "臨床建議與運用..."
    }}
    DO NOT use Markdown. Return ONLY the JSON string.
    """
    
    payload = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            txt = response.json()['candidates'][0]['content']['parts'][0]['text']
            txt = clean_json_text(txt)
            try:
                data = json.loads(txt)
                html_output = f"""
                <div style="font-family: sans-serif; line-height: 1.6; color: #333; background-color: #fff; padding: 10px; border-radius: 8px; border: 1px solid #eee;">
                    <div style="margin-bottom: 15px;">
                        <h4 style="color:#2e86c1; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 5px;">1. 🧪 研究方法 (Methods)</h4>
                        <div style="font-size: 0.95em;">{data.get('methods', '無資料')}</div>
                    </div>
                    <div style="margin-bottom: 15px;">
                        <h4 style="color:#2e86c1; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 5px;">2. 💡 發想緣起 (Rationale)</h4>
                        <div style="font-size: 0.95em;">{data.get('rationale', '無資料')}</div>
                    </div>
                    <div style="margin-bottom: 15px;">
                        <h4 style="color:#2e86c1; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 5px;">3. 📊 結果數據 (Results)</h4>
                        <div style="font-size: 0.95em;">{data.get('results', '無資料')}</div>
                    </div>
                    <div>
                        <h4 style="color:#d35400; margin:0 0 5px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 5px;">4. 🏥 臨床運用 (Implication)</h4>
                        <div style="font-size: 0.95em;">{data.get('implication', '無資料')}</div>
                    </div>
                </div>
                """
                return html_output
            except json.JSONDecodeError: return f"<div style='color:red'>JSON 解析失敗，請重試</div>"
        else: return f"<div style='color:red'>API Error: {response.status_code}</div>"
    except Exception as e: return f"<div style='color:red'>System Error: {str(e)}</div>"

def send_bulk_email(to, pwd, queue):
    if not queue: return False, "清單為空"
    msg = MIMEMultipart()
    msg['From'] = to
    msg['To'] = to
    msg['Subject'] = f"GynOnc 文獻彙報 ({len(queue)}篇) - {datetime.now().strftime('%Y-%m-%d')}"
    
    body = """
    <html><body style="font-family: Arial, sans-serif; color: #333;">
    <h2 style="color: #2c3e50;">🧬 GynOnc 文獻分析報告</h2>
    <p>以下是您精選的文獻深度分析：</p>
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

st.title("💎 GynOnc 文獻系統 v6.2")
st.caption("修復版：解決變數 NameError，確保搜尋功能正常")

# 1. 執行搜尋
if st.session_state.search_trigger:
    with st.status("🔍 正在執行搜尋...", expanded=True) as status:
        q = build_query(final_keywords, final_journals, date_range_query)
        st.write(f"語法: `{q[:50]}...`")
        
        # 這裡的 user_email 現在一定有定義了
        raw_articles = fetch_headers(q, date_params, max_results, user_email)
        
        if raw_articles:
            st.write(f"✅ 找到 {len(raw_articles)} 篇，翻譯標題中...")
            translated_articles = batch_translate_titles_robust(raw_articles, api_key, selected_model_name)
            st.session_state.articles_data = translated_articles
            status.update(label="搜尋完成！", state="complete")
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
                st.markdown(f"<span style='color:#2e86c1; font-size:1.1em;'>{art.get('title_zh', '翻譯中...')}</span>", unsafe_allow_html=True)
                st.caption(f"📖 {art['journal']} | [原文連結]({art['link']})")
            
            with col2:
                btn_key = f"analyze_{art['id']}_{i}"
                if st.button("🔍 詳細分析", key=btn_key):
                    with st.spinner("AI 正在深度閱讀並生成報告..."):
                        if art['id'] not in st.session_state.analysis_cache:
                            report_html = run_deep_analysis_json(art, api_key, selected_model_name)
                            st.session_state.analysis_cache[art['id']] = report_html
                            
                            email_item = {
                                "title": art['title'],
                                "html": f"""
                                <div style="background-color: #f9f9f9; padding: 20px; border-radius: 5px; margin-bottom: 20px;">
                                    <h3 style="margin-top: 0; color: #1a5276;"><a href='{art['link']}' style="text-decoration: none; color: #1a5276;">{art['title']}</a></h3>
                                    <h4 style="margin-top: 5px; color: #2e86c1;">{art.get('title_zh', '')}</h4>
                                    <p style="color: #666; font-size: 0.9em;">📖 {art['journal']}</p>
                                    {report_html}
                                </div>
                                """
                            }
                            if not any(d['title'] == art['title'] for d in st.session_state.email_queue):
                                st.session_state.email_queue.append(email_item)
                                st.rerun()

            if art['id'] in st.session_state.analysis_cache:
                with st.expander("🩺 深度分析報告", expanded=True):
                    st.markdown(st.session_state.analysis_cache[art['id']], unsafe_allow_html=True)
            st.markdown("---")

if getattr(st.session_state, 'trigger_email', False):
    ok, msg = send_bulk_email(user_email, email_password, st.session_state.email_queue)
    if ok:
        st.sidebar.success("✅ 郵件已成功寄出！")
        st.session_state.email_queue = []
        st.session_state.trigger_email = False
        time.sleep(2)
        st.rerun()
    else:
        st.sidebar.error(f"❌ 寄送失敗: {msg}")
        st.session_state.trigger_email = False
