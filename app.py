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
    당신은 유전체 분석가입니다. 제공된 텍스트에서 유전자와 질환 관계를 추출하여 JSON으로만 응답하세요.
    반드시 'nodes'와 'links' 키를 포함하고, 각 노드에는 'id', 'label', 'type', 'desc' 필드를 포함하세요.
    
    형식 예시:
    {{
      "nodes": [{{ "id": "G1", "label": "BRCA1", "type": "gene", "desc": "유방암 위험 증가와 관련된 유전자" }}],
      "links": [{{ "source": "G1", "target": "D1" }}]
    }}
    
    텍스트: {text[:15000]}
    """
    try:
        response = model.generate_content(prompt)
        # JSON 부분만 추출
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group())
        return None
    except Exception as e:
        st.error(f"AI 분석 중 오류 발생: {e}")
        return None

# --- 메인 UI: 다중 파일 업로드 ---
uploaded_files = st.file_uploader("PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

# 업로드된 파일들의 이름 리스트를 추적하여 변경 감지
current_file_names = [f.name for f in uploaded_files] if uploaded_files else []

if uploaded_files and api_key:
    # 파일 구성이 바뀌면(새 파일 추가/삭제) 세션 강제 초기화 및 재분석
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
            # --- KeyError 방지를 위한 안전한 노드/간선 생성 ---
            nodes = []
            raw_nodes = st.session_state.graph_data.get('nodes', [])
            for n in raw_nodes:
                if 'id' in n:
                    # label이 없으면 id를 대신 사용, desc가 없으면 기본값 사용
                    l_text = n.get('label', n['id'])
                    n_type = n.get('type', 'gene')
                    n_color = '#4285F4' if n_type == 'gene' else '#EA4335'
                    nodes.append(Node(id=n['id'], label=l_text, size=25, color=n_color))
            
            edges = []
            raw_links = st.session_state.graph_data.get('links', [])
            for l in raw_links:
                if 'source' in l and 'target' in l:
                    edges.append(Edge(source=l['source'], target=l['target']))

            # 노드가 존재할 때만 그래프 렌더링
            if nodes:
                config = Config(
                    width=900, 
                    height=600, 
                    directed=True, 
                    physics=True, 
                    hierarchical=False,
                    fit_view=True  # 그래프를 화면 중앙에 자동으로 맞춤
                )
                selected_id = agraph(nodes=nodes, edges=edges, config=config)
            else:
                st.warning("그래프를 생성할 유효한 데이터가 없습니다.")
                selected_id = None

        with col2:
            st.markdown("### 🔍 상세 정보")
            if selected_id:
                # 선택된 노드의 상세 정보 찾기
                node_detail = next((n for n in st.session_state.graph_data.get('nodes', []) if str(n.get('id')) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail.get('label', selected_id)}")
                    st.info(f"**설명:** {node_detail.get('desc', '설명 정보가 없습니다.')}")
                else:
                    st.write("상세 정보를 불러올 수 없습니다.")
            else:
                st.write("💡 그래프의 노드를 클릭하면 상세 정보가 여기에 표시됩니다.")

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
            with st.spinner("생각 중..."):
                try:
                    # 문서 내용과 함께 질문 던지기
                    response = model.generate_content(f"문서 통합본 내용: {st.session_state.full_text[:8000]}\n\n질문: {prompt}")
                    st.markdown(response.text)
                    st.session_state.messages.append({"role": "assistant", "content": response.text})
                except Exception as e:
                    st.error(f"채팅 응답 중 오류가 발생했습니다: {e}")