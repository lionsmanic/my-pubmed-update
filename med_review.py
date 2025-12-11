# --- 側邊欄 (優先顯示購物車) ---
with st.sidebar:
    st.header("💎 設定與購物車")
    
    # --- 1. 購物車/寄信區 ---
    if st.session_state.email_queue:
        with st.expander(f"🛒 待寄出清單 ({len(st.session_state.email_queue)}篇)", expanded=True):
            for item in st.session_state.email_queue:
                st.text(f"• {item['title'][:20]}...")
            
            if 'EMAIL_ADDRESS' in st.secrets: user_email = st.secrets['EMAIL_ADDRESS']
            else: user_email = st.text_input("您的 Email", "lionsmanic@gmail.com")
            
            if 'EMAIL_PASSWORD' in st.secrets: email_password = st.secrets['EMAIL_PASSWORD']
            else: email_password = st.text_input("Gmail App Password", type="password")

            if st.button("📩 立即彙整寄出", type="primary"):
                if not email_password:
                    st.error("請輸入密碼")
                else:
                    st.session_state.trigger_email = True
    else:
        st.info("尚無選定文章。請在右側點擊「詳細分析」加入清單。")
    
    st.divider()

    # --- 2. API Key 與模型 (關鍵修正區) ---
    if 'GEMINI_API_KEY' in st.secrets:
        api_key = st.secrets['GEMINI_API_KEY']
        st.success("🔑 API Key 已載入 (Secrets)")
    else:
        api_key = st.text_input("Gemini API Key", type="password")

    # [修正] 預設給一個模型，防止按鈕卡死
    selected_model_name = "gemini-1.5-flash" 
    
    if api_key:
        with st.spinner("正在連線 Google 取得模型清單..."):
            available_models = get_available_models(api_key)
        
        if available_models:
            st.success(f"✅ 連線成功！偵測到 {len(available_models)} 個模型")
            default_ix = 0
            if 'gemini-1.5-flash' in available_models: default_ix = available_models.index('gemini-1.5-flash')
            elif 'gemini-pro' in available_models: default_ix = available_models.index('gemini-pro')
            selected_model_name = st.selectbox("選擇模型:", available_models, index=default_ix)
        else:
            # 即使抓不到，也顯示紅字但不鎖死按鈕
            st.warning("⚠️ 無法自動取得模型清單 (可能是 Key 權限或網路問題)。\n已強制使用預設值：gemini-1.5-flash")
            selected_model_name = st.text_input("手動輸入模型名稱", "gemini-1.5-flash")

    st.divider()
    
    # --- 3. 搜尋條件 ---
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

    # --- 4. 時間與數量 ---
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
    
    # [修正] 移除 disabled 限制，讓您隨時可按
    if st.button("🚀 極速搜尋", type="primary"):
        if not api_key:
            st.error("❌ 請先輸入 API Key 才能搜尋！")
        else:
            st.session_state.articles_data = []
            st.session_state.analysis_cache = {}
            st.session_state.search_trigger = True
