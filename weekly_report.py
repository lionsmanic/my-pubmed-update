import os
import smtplib
import google.generativeai as genai
from Bio import Entrez
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, timedelta

# --- 設定環境變數 (由 GitHub Secrets 提供) ---
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
EMAIL_ADDRESS = os.environ.get("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.environ.get("EMAIL_PASSWORD") # 這是 Gmail 應用程式密碼
Entrez.email = EMAIL_ADDRESS

# --- 設定搜尋參數 ---
# 搜尋過去 7 天的文章
SEARCH_QUERY = '("Ovarian Neoplasms"[Mesh] OR "Uterine Cervical Neoplasms"[Mesh] OR "Gynecologic Neoplasms"[Mesh])'
MAX_RESULTS = 5

def fetch_recent_articles():
    print("正在搜尋 PubMed...")
    # reldate=7 代表搜尋最近 7 天
    handle = Entrez.esearch(db="pubmed", term=SEARCH_QUERY, reldate=7, retmax=MAX_RESULTS, sort="date")
    record = Entrez.read(handle)
    id_list = record["IdList"]
    
    if not id_list:
        return []
        
    handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
    articles = Entrez.read(handle)
    
    parsed_data = []
    for art in articles['PubmedArticle']:
        try:
            citation = art['MedlineCitation']
            title = citation['Article']['ArticleTitle']
            journal = citation['Article']['Journal']['Title']
            
            if 'Abstract' in citation['Article']:
                abstract = " ".join([str(x) for x in citation['Article']['Abstract']['AbstractText']])
            else:
                abstract = "No Abstract"
                
            ids = art['PubmedData']['ArticleIdList']
            doi = next((item for item in ids if item.attributes['IdType'] == 'doi'), None)
            link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{ids[0]}/"
            
            parsed_data.append({"title": title, "abstract": abstract, "link": link, "journal": journal})
        except:
            continue
    return parsed_data

def analyze_and_format_html(articles):
    genai.configure(api_key=GEMINI_API_KEY)
    model = genai.GenerativeModel('gemini-1.5-flash')
    
    email_body = f"<h2>🧬 婦科腫瘤每週文獻快報 ({datetime.now().strftime('%Y-%m-%d')})</h2><hr>"
    
    for i, art in enumerate(articles):
        print(f"正在分析第 {i+1} 篇: {art['title'][:20]}...")
        prompt = f"""
        請閱讀以下摘要，並用繁體中文寫一段簡短的分析（HTML 格式）。
        重點放在：1.背景與結果 2.臨床意義(Clinical Implication) 3.對婦癌醫師的啟發。
        請直接輸出 HTML 標籤（如 <p>, <b>, <ul>），不要用 Markdown。
        
        標題: {art['title']}
        期刊: {art['journal']}
        摘要: {art['abstract']}
        """
        try:
            response = model.generate_content(prompt)
            analysis_html = response.text.replace("```html", "").replace("```", "") # 清理可能的多餘標籤
            
            # 組合單篇文章的 HTML 區塊
            email_body += f"""
            <div style="margin-bottom: 30px; padding: 15px; background-color: #f9f9f9; border-left: 5px solid #2e86c1;">
                <h3 style="color: #1a5276; margin-top: 0;"><a href="{art['link']}">{art['title']}</a></h3>
                <p style="font-size: 0.9em; color: #666;">📖 {art['journal']}</p>
                <div>{analysis_html}</div>
            </div>
            """
        except Exception as e:
            print(f"分析失敗: {e}")
            
    email_body += "<p style='color: #888; font-size: 0.8em;'>本郵件由 Python 自動排程發送，內容由 Gemini AI 生成供參考。</p>"
    return email_body

def send_email(content_html):
    msg = MIMEMultipart()
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = EMAIL_ADDRESS
    msg['Subject'] = f"週一晨報：最新的婦科腫瘤研究 ({datetime.now().strftime('%m/%d')})"
    
    msg.attach(MIMEText(content_html, 'html'))
    
    try:
        # 使用 Gmail SMTP
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
        print("✅ 郵件發送成功！")
    except Exception as e:
        print(f"❌ 郵件發送失敗: {e}")

if __name__ == "__main__":
    articles = fetch_recent_articles()
    if articles:
        print(f"找到 {len(articles)} 篇文章，開始處理...")
        html_content = analyze_and_format_html(articles)
        send_email(html_content)
    else:
        print("本週沒有符合條件的新文章。")
