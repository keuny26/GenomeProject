import streamlit as st
import fitz  # PyMuPDF
import json
import re
import google.generativeai as genai
from streamlit_agraph import agraph, Node, Edge, Config

# --- 페이지 설정 ---
st.set_page_config(page_title="GenomeGraph AI", layout="wide")
st.title("🧬 GenomeGraph AI (Smart Integration)")

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

# --- 분석 함수 ---
def analyze_graph_with_ai(text):
    if not model: return None
    prompt = f"""
    유전체 분석가로서 텍스트에서 유전자와 질환 관계를 추출해 JSON으로만 답하세요. 
    반드시 'nodes'와 'links'를 포함하고 각 노드에 'desc'를 넣으세요.
    텍스트: {text[:15000]}
    """
    try:
        response = model.generate_content(prompt)
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        return json.loads(json_match.group()) if json_match else None
    except: return None

# --- 메인 UI: 다중 파일 업로드 ---
uploaded_files = st.file_uploader("PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

# 핵심 변경 사항: 업로드된 파일들의 이름 리스트를 추적하여 변경 감지
current_file_names = [f.name for f in uploaded_files] if uploaded_files else []

if uploaded_files and api_key:
    # 파일 구성이 바뀌면(새 파일 추가/삭제) 세션 강제 초기화
    if "last_files" not in st.session_state or st.session_state.last_files != current_file_names:
        with st.spinner("새로운 파일을 포함하여 통합 분석 중..."):
            combined_text = ""
            for uploaded_file in uploaded_files:
                doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                combined_text += f"\n[Document: {uploaded_file.name}]\n"
                combined_text += " ".join([page.get_text() for page in doc])
            
            st.session_state.full_text = combined_text
            st.session_state.last_files = current_file_names  # 파일 리스트 업데이트
            st.session_state.graph_data = analyze_graph_with_ai(st.session_state.full_text)
            st.session_state.messages = []

    # 2. 그래프 영역
    if st.session_state.get("graph_data"):
        st.subheader("🧬 통합 지식 그래프")
        
        col1, col2 = st.columns([3, 1])
        
        with col1:
            nodes = [Node(id=n['id'], label=n['label'], size=25, color=('#4285F4' if n.get('type') == 'gene' else '#EA4335')) 
                     for n in st.session_state.graph_data.get('nodes', [])]
            edges = [Edge(source=l['source'], target=l['target']) for l in st.session_state.graph_data.get('links', [])]
            
            # 그래프 설정 개선 (중앙 정렬 및 물리 엔진 최적화)
            config = Config(
                width=900, 
                height=600, 
                directed=True, 
                physics=True, 
                hierarchical=False, # 유전체 관계는 계층보다 네트워크 형태가 적합
                fit_view=True       # 그래프를 화면 중앙에 자동으로 맞춤
            )
            selected_id = agraph(nodes=nodes, edges=edges, config=config)

        with col2:
            st.markdown("### 🔍 상세 정보")
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data['nodes'] if n['id'] == selected_id), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail['label']}")
                    st.info(f"**설명:** {node_detail.get('desc', '설명 없음')}")
            else:
                st.write("💡 노드를 클릭하세요.")

    st.divider()

    # 3. 채팅 영역
    st.subheader("💬 통합 분석 채팅")
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("질문하세요."):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)
        with st.chat_message("assistant"):
            response = model.generate_content(f"문서 통합본: {st.session_state.full_text[:8000]}\n\n질문: {prompt}")
            st.markdown(response.text)
            st.session_state.messages.append({"role": "assistant", "content": response.text})