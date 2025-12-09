import streamlit as st
import google.generativeai as genai
from Bio import Entrez
import datetime

# --- 頁面設定 ---
st.set_page_config(
    page_title="婦科腫瘤文獻智慧分析",
    page_icon="🧬",
    layout="wide"
)

# 修改側邊欄程式碼片段
with st.sidebar:
    st.header("⚙️ 設定控制台")
    
    # 嘗試從 Secrets 讀取，如果沒有才讓使用者輸入
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("✅ 已從系統讀取 API Key")
    else:
        api_key = st.text_input("請輸入 Google Gemini API Key", type="password")

    if 'EMAIL_ADDRESS' in st.secrets:
        email_input = st.secrets['EMAIL_ADDRESS']
    else:
        email_input = st.text_input("Email", "lionsmanic@gmail.com")
    
    st.divider()
    
    # 2. 搜尋參數
    st.subheader("🔍 搜尋條件")
    # 預設一些婦癌關鍵字
    default_query = '("Ovarian Neoplasms"[Mesh] OR "Uterine Cervical Neoplasms"[Mesh]) AND "2024"[Date - Publication]'
    query = st.text_area("PubMed 搜尋語法 (支援布林邏輯)", value=default_query, height=100)
    
    st.info("💡 提示：您可以輸入 'Lancet Oncol[Journal]' 來鎖定特定期刊。")
    
    max_results = st.slider("分析篇數 (建議 3-5 篇以節省時間)", 1, 10, 3)
    
    # 按鈕
    start_btn = st.button("🚀 開始搜尋與分析", type="primary")

# --- 核心函數：抓取 PubMed ---
def fetch_pubmed_articles(query, max_results, email):
    Entrez.email = email
    try:
        # 1. 搜尋 ID
        search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="date")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        if not id_list:
            return []

        # 2. 抓取詳細內容
        fetch_handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        articles_data = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        parsed_articles = []
        for article in articles_data['PubmedArticle']:
            try:
                citation = article['MedlineCitation']
                title = citation['Article']['ArticleTitle']
                
                # 處理摘要列表
                if 'Abstract' in citation['Article']:
                    abstract_parts = citation['Article']['Abstract']['AbstractText']
                    abstract = " ".join([str(part) for part in abstract_parts])
                else:
                    abstract = "無摘要 (No Abstract Available)"
                
                # 抓取期刊與年份
                journal = citation['Article']['Journal']['Title']
                pub_date = citation['Article']['Journal']['JournalIssue']['PubDate']
                date_str = f"{pub_date.get('Year', '')} {pub_date.get('Month', '')}"
                
                # 抓取 DOI 連結
                ids = article['PubmedData']['ArticleIdList']
                doi = next((item for item in ids if item.attributes['IdType'] == 'doi'), None)
                link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{ids[0]}/"

                parsed_articles.append({
                    "title": title,
                    "abstract": abstract,
                    "journal": journal,
                    "date": date_str,
                    "link": link
                })
            except Exception as e:
                continue # 跳過格式錯誤的文章
                
        return parsed_articles

    except Exception as e:
        st.error(f"PubMed 連線錯誤: {e}")
        return []

# --- 核心函數：Gemini 分析 ---
def analyze_article(article, api_key):
    genai.configure(api_key=api_key)
    # 使用 Flash 模型速度較快且便宜，若需要更深度的推理可改用 pro
    model = genai.GenerativeModel('gemini-1.5-flash') 
    
    prompt = f"""
    你現在是一位權威的「婦科腫瘤學教授」與臨床醫師。請閱讀以下這篇醫學文獻的摘要，並為你的主治醫師團隊用「繁體中文」做重點解讀。
    
    【文獻資訊】
    標題: {article['title']}
    期刊: {article['journal']}
    摘要: {article['abstract']}
    
    【請依序輸出以下區塊，並使用 Markdown 格式排版】：

    ### 1. 📝 文獻快報 (Structured Summary)
    請簡明扼要地整理：
    * **Background (背景)**: 
    * **Methods (方法)**: 
    * **Results (主要結果)**: (請包含重要的統計數據，如 P值、HR、OR等)
    * **Conclusion (結論)**: 

    ### 2. 💡 發想緣起 (Origin of the Idea)
    (請根據背景推論：為什麼作者想做這個研究？是為了解決什麼過去臨床上的痛點、爭議或是補足哪塊證據？)

    ### 3. 🏥 臨床可行運用 (Clinical Application)
    (這對我們目前的臨床實踐有什麼直接影響？是否支持改變現有的治療策略？例如手術方式、化療藥物選擇或篩檢流程？若尚不可行，請說明原因。)

    ### 4. 🚀 婦癌醫師的研究機遇 (Future Directions for GynOnc)
    (針對婦科腫瘤醫師，這篇研究啟發了什麼後續方向？有沒有我們可以在本地醫院利用現有病歷資料進行驗證的題目？或是延伸的子題？)
    """
    
    try:
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        return f"AI 分析失敗: {e}"

# --- 主畫面邏輯 ---
st.title("🧬 GynOnc 醫學文獻智慧分析助手")
st.markdown("專為婦科腫瘤醫師設計，自動抓取 PubMed 並生成臨床應用導向的分析報告。")

if start_btn:
    if not api_key:
        st.warning("⚠️ 請先在側邊欄輸入 Gemini API Key")
    else:
        with st.status("🔄 正在執行任務中...", expanded=True) as status:
            
            # 1. 搜尋
            st.write("📡 連接 PubMed 資料庫搜尋中...")
            articles = fetch_pubmed_articles(query, max_results, email_input)
            
            if not articles:
                status.update(label="❌ 搜尋不到結果，請檢查關鍵字", state="error")
            else:
                st.write(f"✅ 成功抓取 {len(articles)} 篇文章，開始 AI 閱讀分析...")
                
                # 建立一個容器來放結果
                results_container = st.container()
                
                for i, article in enumerate(articles):
                    st.write(f"🤖 正在分析第 {i+1} 篇: {article['title'][:30]}...")
                    
                    # 呼叫 Gemini
                    analysis = analyze_article(article, api_key)
                    
                    # 顯示結果
                    with results_container:
                        st.markdown("---")
                        st.subheader(f"#{i+1} {article['title']}")
                        st.caption(f"📖 {article['journal']} | 🗓️ {article['date']}")
                        st.markdown(f"🔗 [點擊閱讀原文]({article['link']})")
                        
                        with st.expander("查看原始摘要 (English Abstract)"):
                            st.text(article['abstract'])
                        
                        # AI 輸出區塊 - 重點樣式
                        st.info("🤖 **Gemini 教授的分析報告**")
                        st.markdown(analysis)
                
                status.update(label="🎉 所有文獻分析完成！", state="complete")
