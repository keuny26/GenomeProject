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

# --- 분석 함수 (멀티 파일 통합 최적화) ---
def analyze_graph_with_ai(text):
    if not model: return None
    
    # 2개 이상의 문서를 모두 읽도록 하는 강력한 지시사항
    prompt = f"""
    당신은 전문 유전체 분석가입니다. 
    제공된 텍스트는 여러 개의 유전체 보고서가 통합된 내용입니다.

    [작업 목표]
    1. **모든 문서 포함**: 제공된 모든 텍스트([Document: 파일명] 태그로 구분됨)를 훑어보고 각 문서의 유전자-질환 관계를 추출하세요.
    2. **통합 및 연결**: 여러 문서에서 공통으로 등장하는 노드는 하나로 합치고, 각 문서의 분석 내용을 'desc' 필드에 종합하세요.
    3. **JSON 출력**: 반드시 아래 구조를 지킨 JSON 데이터만 응답하세요.

    {{
      "nodes": [
        {{ "id": "unique_id", "label": "유전자/질환명", "type": "gene 또는 disease", "desc": "상세 내용 요약" }}
      ],
      "links": [
        {{ "source": "node_id", "target": "node_id" }}
      ]
    }}

    텍스트 (최대 20,000자):
    {text[:20000]}
    """
    
    try:
        response = model.generate_content(prompt)
        # JSON 블록만 추출하는 정규표현식
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group())
        return None
    except Exception as e:
        if "429" in str(e):
            st.error("API 할당량을 초과했습니다. 약 1분 후 다시 시도하세요.")
        else:
            st.error(f"AI 분석 중 오류 발생: {e}")
        return None

# --- 메인 UI: 다중 파일 업로드 ---
uploaded_files = st.file_uploader("PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

if uploaded_files and api_key:
    # 쿼터 보호를 위한 분석 시작 버튼
    if st.button("🧬 통합 분석 시작 (API 호출)"):
        with st.spinner("여러 문서의 데이터를 통합 분석 중입니다..."):
            combined_text = ""
            for uploaded_file in uploaded_files:
                doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                # 각 문서의 시작임을 AI에게 알림
                combined_text += f"\n\n[Document: {uploaded_file.name}]\n"
                combined_text += " ".join([page.get_text() for page in doc])
            
            st.session_state.full_text = combined_text
            st.session_state.last_files = [f.name for f in uploaded_files]
            # AI 분석 호출
            st.session_state.graph_data = analyze_graph_with_ai(st.session_state.full_text)
            st.session_state.messages = []
            
            if st.session_state.graph_data:
                st.success(f"성공: {len(uploaded_files)}개의 파일 분석 완료!")
            else:
                st.error("AI가 데이터를 추출하지 못했습니다. 파일 내용을 확인해주세요.")

    # 2. 그래프 영역
    if st.session_state.get("graph_data"):
        st.subheader("🧬 통합 지식 그래프")
        st.info("💡 줌인/아웃이 가능하며, 노드를 클릭하면 상세 정보가 표시됩니다.")
        
        col1, col2 = st.columns([3, 1])
        
        with col1:
            nodes = []
            raw_nodes = st.session_state.graph_data.get('nodes', [])
            for n in raw_nodes:
                if 'id' in n:
                    l_text = n.get('label', n['id'])
                    n_type = n.get('type', '').lower()
                    # 유전자는 파란색, 질환은 빨간색으로 구분
                    n_color = '#4285F4' if 'gene' in n_type else '#EA4335'
                    nodes.append(Node(id=n['id'], label=l_text, size=25, color=n_color))
            
            edges = []
            raw_links = st.session_state.graph_data.get('links', [])
            for l in raw_links:
                if 'source' in l and 'target' in l:
                    edges.append(Edge(source=l['source'], target=l['target']))

            if nodes:
                # fit_view=True로 설정하여 항상 그래프가 중앙에 오도록 함
                config = Config(
                    width=900, 
                    height=600, 
                    directed=True, 
                    physics=True, 
                    fit_view=True,
                    nodeHighlightBehavior=True,
                    highlightColor="#F79767"
                )
                selected_id = agraph(nodes=nodes, edges=edges, config=config)
            else:
                st.warning("표시할 노드가 없습니다.")
                selected_id = None

        with col2:
            st.markdown("### 🔍 상세 정보")
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data.get('nodes', []) if str(n.get('id')) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail.get('label', selected_id)}")
                    st.markdown(f"**설명:**\n{node_detail.get('desc', '설명 정보가 없습니다.')}")
            else:
                st.write("💡 그래프에서 노드를 클릭하세요.")

        st.divider()

        # 3. 채팅 영역
        st.subheader("💬 통합 분석 채팅")
        for message in st.session_state.messages:
            with st.chat_message(message["role"]):
                st.markdown(message["content"])

        if prompt := st.chat_input("이 문서들의 공통점이나 특정 유전자에 대해 물어보세요."):
            st.session_state.messages.append({"role": "user", "content": prompt})
            with st.chat_message("user"):
                st.markdown(prompt)
            with st.chat_message("assistant"):
                with st.spinner("통합 문서 확인 중..."):
                    try:
                        # 채팅 시에도 통합 텍스트 기반으로 답변
                        response = model.generate_content(f"당신은 유전체 전문가입니다. 다음 통합 보고서 내용을 바탕으로 답변하세요.\n내용: {st.session_state.full_text[:8000]}\n\n질문: {prompt}")
                        st.markdown(response.text)
                        st.session_state.messages.append({"role": "assistant", "content": response.text})
                    except Exception as e:
                        st.error(f"채팅 중 오류 발생: {e}")