import streamlit as st
import fitz  # PyMuPDF
import json
import re
import google.generativeai as genai
from streamlit_agraph import agraph, Node, Edge, Config

# --- 페이지 설정 ---
st.set_page_config(page_title="GenomeGraph AI", layout="wide")
st.title("🧬 GenomeGraph AI (Streamlit)")

# --- API 키 설정 ---
if "GEMINI_API_KEY" in st.secrets:
    api_key = st.secrets["GEMINI_API_KEY"]
else:
    st.sidebar.title("설정")
    api_key = st.sidebar.text_input("Gemini API Key를 입력하세요", type="password")

# --- 모델 초기화 ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        available_models = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]
        
        target_model = 'gemini-1.5-flash'
        if f'models/{target_model}' in available_models:
            model = genai.GenerativeModel(target_model)
        elif available_models:
            model = genai.GenerativeModel(available_models[0].replace('models/', ''))
        
        if model:
            st.sidebar.success(f"연결됨: {model.model_name}")
    except Exception as e:
        st.error(f"API 설정 오류: {e}")

# --- 분석 함수 (그래프 데이터 생성) ---
def analyze_graph_with_ai(text):
    if not model: return None
    prompt = f"""
    당신은 유전체 데이터 분석가입니다. 아래 텍스트에서 유전자와 질환 관계를 추출하여 JSON으로만 응답하세요.
    각 노드에는 반드시 'desc'(상세 설명) 필드를 포함시켜주세요.

    형식:
    {{
      "nodes": [{{ "id": "ID", "label": "이름", "type": "gene/disease", "desc": "상세 분석 내용" }}],
      "links": [{{ "source": "ID", "target": "ID" }}]
    }}
    텍스트: {text[:10000]}
    """
    try:
        response = model.generate_content(prompt)
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        return json.loads(json_match.group()) if json_match else None
    except: return None

# --- 메인 UI ---
uploaded_file = st.file_uploader("PDF 보고서를 업로드하세요", type="pdf")

if uploaded_file and api_key:
    # 1. 데이터 로드 (세션 상태 저장)
    if "full_text" not in st.session_state:
        with st.spinner("PDF 분석 중..."):
            doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
            st.session_state.full_text = " ".join([page.get_text() for page in doc])
            st.session_state.graph_data = analyze_graph_with_ai(st.session_state.full_text)
            st.session_state.messages = [] # 채팅 초기화

    # 2. 그래프 및 상세 정보 레이아웃
    if st.session_state.get("graph_data"):
        st.subheader("🧬 유전체 지식 그래프 및 상세 정보")
        
        col1, col2 = st.columns([3, 1])
        
        with col1:
            nodes = [Node(id=n['id'], label=n['label'], size=20, color=('#4285F4' if n.get('type') == 'gene' else '#EA4335')) 
                     for n in st.session_state.graph_data.get('nodes', [])]
            edges = [Edge(source=l['source'], target=l['target']) for l in st.session_state.graph_data.get('links', [])]
            
            config = Config(width=800, height=500, directed=True, physics=True)
            # 노드 클릭 시 selected_id에 해당 노드의 ID가 저장됨
            selected_id = agraph(nodes=nodes, edges=edges, config=config)

        with col2:
            st.markdown("### 🔍 상세 정보")
            if selected_id:
                # 선택된 노드 찾기
                node_detail = next((n for n in st.session_state.graph_data['nodes'] if n['id'] == selected_id), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail['label']}")
                    st.info(f"**설명:** {node_detail.get('desc', '이 노드에 대한 상세 설명이 없습니다.')}")
                else:
                    st.warning("데이터 동기화 중... 잠시만 기다려주세요.")
            else:
                st.write("💡 그래프의 동그라미(노드)를 클릭하면 상세 분석 내용이 여기에 표시됩니다.")

    st.divider()

    # 3. AI 채팅창 영역
    st.subheader("💬 AI 분석가와 대화하기")
    
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("보고서 내용에 대해 질문하세요."):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

        with st.chat_message("assistant"):
            with st.spinner("답변 생성 중..."):
                response = model.generate_content(f"문서 내용: {st.session_state.full_text[:8000]}\n\n질문: {prompt}")
                st.markdown(response.text)
                st.session_state.messages.append({"role": "assistant", "content": response.text})