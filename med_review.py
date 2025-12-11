def run_deep_analysis(art, key, model_name):
    """單篇深度分析 (修正版：強制 HTML 格式)"""
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={key}"
    headers = {'Content-Type': 'application/json'}
    
    # 修正後的 Prompt：強制要求統一的 HTML 結構
    prompt_text = f"""
    角色：資深婦癌權威醫師。
    任務：針對以下這篇論文進行詳細的學術點評。
    
    標題：{art['title']}
    摘要：{art['abstract']}
    
    【輸出格式要求】：
    1. 請輸出 **純 HTML 代碼**。
    2. **嚴禁**使用 Markdown (如 ##, **, - 等符號)。
    3. 所有標題都必須使用 <h4 style="color:#2e86c1;"> 標籤包覆。
    4. 所有內容必須包在一個最外層的 <div> 標籤內。
    
    請依照以下模板輸出：
    <div style="font-family: sans-serif; line-height: 1.6;">
        <h4 style="color:#2e86c1; margin-top:0; border-bottom: 2px solid #eee; padding-bottom: 5px;">1. 🧪 研究方法 (Methods)</h4>
        <p>請簡述 Study Design, Patient Population, Intervention。</p>
        
        <h4 style="color:#2e86c1; border-bottom: 2px solid #eee; padding-bottom: 5px;">2. 💡 發想緣起 (Rationale)</h4>
        <p>推測作者為何進行此研究？解決了什麼臨床痛點？</p>
        
        <h4 style="color:#2e86c1; border-bottom: 2px solid #eee; padding-bottom: 5px;">3. 📊 結果數據 (Results)</h4>
        <ul>
            <li>關鍵 P-value: ...</li>
            <li>HR / OR: ...</li>
        </ul>
        
        <h4 style="color:#d35400; border-bottom: 2px solid #eee; padding-bottom: 5px;">4. 🏥 臨床運用與結論 (Clinical Implication)</h4>
        <p>這對婦癌臨床實踐有何具體改變或建議？</p>
    </div>
    """
    
    payload = {"contents": [{"parts": [{"text": prompt_text}]}]}
    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            txt = response.json()['candidates'][0]['content']['parts'][0]['text']
            # 清理可能的多餘標記
            clean_html = txt.replace("```html", "").replace("```", "").strip()
            return clean_html
        else: 
            return f"<div style='color:red'>分析失敗: {response.text}</div>"
    except Exception as e: 
        return f"<div style='color:red'>錯誤: {str(e)}</div>"
