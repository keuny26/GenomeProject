import streamlit as st
import fitz  # PyMuPDF
import json
import re
import google.generativeai as genai
from streamlit_agraph import agraph, Node, Edge, Config

# --- 페이지 설정 ---
st.set_page_config(page_title="GenomeGraph AI", layout="wide")
st.title("🧬 GenomeGraph AI (Clean Labels)")

# --- 세션 상태 초기화 (AttributeError 방지) ---
if "messages" not in st.session_state:
    st.session_state.messages = []
if "full_text" not in st.session_state:
    st.session_state.full_text = ""
if "graph_data" not in st.session_state:
    st.session_state.graph_data = None

# --- API 키 설정 ---
if "GEMINI_API_KEY" in st.secrets:
    api_key = st.secrets["GEMINI_API_KEY"]
else:
    st.sidebar.title("설정")
    api_key = st.sidebar.text_input("Gemini API Key를 입력하세요", type="password")

# --- 모델 초기화 (자동 감지 로직) ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        available_models = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]
        target_model_name = 'models/gemini-1.5-flash'
        
        if target_model_name in available_models:
            model = genai.GenerativeModel('gemini-1.5-flash')
        elif available_models:
            fallback = available_models[0].replace('models/', '')
            model = genai.GenerativeModel(fallback)
            st.sidebar.warning(f"Flash 모델 미지원으로 {fallback} 모델에 연결되었습니다.")
        
        if model:
            st.sidebar.success(f"연결됨: {model.model_name}")
    except Exception as e:
        st.error(f"API/모델 설정 오류: {e}")

# --- 분석 함수 ---
def analyze_graph_with_ai(text):
    if not model: return None
    prompt = f"""
    당신은 전문 유전체 분석가입니다. 제공된 텍스트에서 유전자와 질환 관계를 추출하세요.
    1. 모든 노드에 'source_file' 필드를 추가하여 출처를 기록하세요.
    2. 공통 노드는 'source_file'을 "Common"으로 하세요.
    3. 반드시 JSON 형식으로만 응답하세요.
    텍스트: {text[:20000]}
    """
    try:
        response = model.generate_content(prompt)
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group())
        return None
    except Exception as e:
        st.error(f"AI 분석 중 오류: {e}")
        return None

# --- 메인 UI: 다중 파일 업로드 ---
uploaded_files = st.file_uploader("PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

if uploaded_files and api_key:
    if st.button("🧬 파일별 통합 분석 시작"):
        with st.spinner("데이터 분석 중..."):
            combined_text = ""
            for uploaded_file in uploaded_files:
                doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                combined_text += f"\n\n[Document: {uploaded_file.name}]\n"
                combined_text += " ".join([page.get_text() for page in doc])
            
            st.session_state.full_text = combined_text
            st.session_state.graph_data = analyze_graph_with_ai(st.session_state.full_text)
            st.session_state.messages = [] 
            st.success("분석이 완료되었습니다!")

    # 2. 그래프 영역
    if st.session_state.graph_data:
        st.subheader("🧬 출처별 통합 지식 그래프")
        col1, col2 = st.columns([3, 1])
        
        file_names = [f.name for f in uploaded_files]
        color_palette = ["#4285F4", "#34A853", "#FBBC05", "#8E44AD", "#F39C12", "#16A085"]
        color_map = {name: color_palette[i % len(color_palette)] for i, name in enumerate(file_names)}
        color_map["Common"] = "#EA4335" 

        with col1:
            nodes = []
            for n in st.session_state.graph_data.get('nodes', []):
                src = n.get('source_file', 'Unknown')
                n_color = color_map.get(src, "#999999")
                n_size = 35 if src == "Common" else 25
                
                # --- ✅ 수정 포인트: [pdf이름] 태그를 삭제하고 순수 이름만 표시 ---
                # 공통 노드일 때만 구분을 위해 별표(⭐)를 유지하고, 아니면 레이블 그대로 사용
                clean_label = f"⭐ {n.get('label')}" if src == "Common" else n.get('label')
                
                nodes.append(Node(id=n['id'], label=clean_label, size=n_size, color=n_color))
            
            edges = [Edge(source=l['source'], target=l['target']) for l in st.session_state.graph_data.get('links', [])]

            if nodes:
                config = Config(width=900, height=600, directed=True, physics=True, fit_view=True, panAndZoom=True)
                selected_id = agraph(nodes=nodes, edges=edges, config=config)
                
                if st.button("🎯 그래프 중앙 정렬"):
                    st.rerun()
            else:
                st.warning("분석 결과가 없습니다.")
                selected_id = None

        with col2:
            st.markdown("### 🎨 범례 및 상세 정보")
            for src, color in color_map.items():
                st.markdown(f"<span style='color:{color}'>●</span> **{src}**", unsafe_allow_html=True)
            
            st.divider()
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data['nodes'] if str(n['id']) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail.get('label')}")
                    st.info(f"**출처:** {node_detail.get('source_file')}")
                    st.write(f"**분석 상세:**\n{node_detail.get('desc', '설명 없음')}")

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
            try:
                res = model.generate_content(f"내용 요약: {st.session_state.full_text[:8000]}\n질문: {prompt}")
                st.markdown(res.text)
                st.session_state.messages.append({"role": "assistant", "content": res.text})
            except Exception as e:
                st.error(f"오류: {e}")