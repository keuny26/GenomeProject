import streamlit as st
import fitz  # PyMuPDF
import json
import re
import google.generativeai as genai
from streamlit_agraph import agraph, Node, Edge, Config

# --- 페이지 설정 ---
st.set_page_config(page_title="GenomeGraph AI", layout="wide")
st.title("🧬 GenomeGraph AI (Smart & Source-Specific)")

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
        model = genai.GenerativeModel('gemini-1.5-flash')
        if model:
            st.sidebar.success(f"연결됨: {model.model_name}")
    except Exception as e:
        st.error(f"API 설정 오류: {e}")

# --- 분석 함수 (출처 태깅 강화) ---
def analyze_graph_with_ai(text):
    if not model: return None
    
    prompt = f"""
    당신은 전문 유전체 분석가입니다. 제공된 텍스트는 여러 보고서의 합본입니다.
    
    [작업 지침]
    1. **출처 명시**: 모든 노드에 'source_file' 필드를 추가하여 어떤 문서([Document: 파일명])에서 추출되었는지 기록하세요.
    2. **공통 노드 처리**: 2개 이상의 문서에서 발견된 유전자/질환은 'source_file' 값을 "Common"으로 지정하세요.
    3. **데이터 구조**: 반드시 아래 JSON 형식을 유지하세요.
    
    {{
      "nodes": [
        {{ "id": "ID", "label": "명칭", "type": "gene/disease", "source_file": "파일명 또는 Common", "desc": "상세 분석" }}
      ],
      "links": [{{ "source": "ID", "target": "ID" }}]
    }}

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
        with st.spinner("파일별 출처를 구분하여 통합 분석 중..."):
            combined_text = ""
            for uploaded_file in uploaded_files:
                doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                combined_text += f"\n\n[Document: {uploaded_file.name}]\n"
                combined_text += " ".join([page.get_text() for page in doc])
            
            st.session_state.full_text = combined_text
            st.session_state.graph_data = analyze_graph_with_ai(st.session_state.full_text)
            st.session_state.messages = []
            st.success("출처별 통합 분석이 완료되었습니다!")

    # 2. 그래프 영역
    if st.session_state.get("graph_data"):
        st.subheader("🧬 출처별 통합 지식 그래프")
        st.info("💡 줌인/아웃이 가능하며, 노드를 클릭하면 상세 정보가 표시됩니다. 노드를 잃어버리면 아래 '🎯 그래프 중앙 정렬'을 누르세요.")
        
        col1, col2 = st.columns([3, 1])
        
        # 색상 팔레트 정의
        color_palette = ["#4285F4", "#34A853", "#FBBC05", "#8E44AD", "#F39C12", "#16A085"]
        file_names = [f.name for f in uploaded_files]
        color_map = {name: color_palette[i % len(color_palette)] for i, name in enumerate(file_names)}
        color_map["Common"] = "#EA4335" 

        with col1:
            nodes = []
            raw_nodes = st.session_state.graph_data.get('nodes', [])
            for n in raw_nodes:
                src = n.get('source_file', 'Unknown')
                n_color = color_map.get(src, "#999999")
                n_size = 35 if src == "Common" else 25
                
                nodes.append(Node(
                    id=n['id'], 
                    label=f"[{src}] {n.get('label')}" if src != "Common" else f"⭐ {n.get('label')}",
                    size=n_size, 
                    color=n_color
                ))
            
            edges = [Edge(source=l['source'], target=l['target']) for l in st.session_state.graph_data.get('links', [])]

            if nodes:
                # --- 화면 이동 및 정렬 설정 최적화 ---
                config = Config(
                    width=900, 
                    height=600, 
                    directed=True, 
                    physics=True, 
                    fit_view=True, # 초기 실행 시 중앙 정렬
                    nodeHighlightBehavior=True,
                    highlightColor="#F79767",
                    panAndZoom=True, # 마우스 드래그 및 휠 지원
                    staticGraph=False
                )
                selected_id = agraph(nodes=nodes, edges=edges, config=config)
                
                # 원이 화면 밖으로 나갔을 때 중앙으로 불러오는 구원 버튼
                if st.button("🎯 그래프 중앙 정렬 (Reset View)"):
                    st.rerun() 
            else:
                st.warning("분석된 데이터가 없습니다.")
                selected_id = None

        with col2:
            st.markdown("### 🎨 범례 및 정보")
            for src, color in color_map.items():
                st.markdown(f"<span style='color:{color}'>●</span> **{src}**", unsafe_allow_html=True)
            
            st.divider()
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data.get('nodes', []) if str(n.get('id')) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail.get('label')}")
                    st.info(f"**출처:** {node_detail.get('source_file')}")
                    st.markdown(f"**분석:** {node_detail.get('desc')}")
            else:
                st.write("💡 노드를 클릭하면 상세 내용이 표시됩니다.")

    st.divider()

    # 3. 채팅 영역
    st.subheader("💬 통합 분석 채팅")
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("이 데이터들에 대해 궁금한 점을 질문하세요."):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)
        with st.chat_message("assistant"):
            with st.spinner("통합 보고서 분석 중..."):
                try:
                    response = model.generate_content(f"다음 통합된 유전체 문서 내용을 바탕으로 답변하세요.\n\n{st.session_state.full_text[:8000]}\n\n질문: {prompt}")
                    st.markdown(response.text)
                    st.session_state.messages.append({"role": "assistant", "content": response.text})
                except Exception as e:
                    st.error(f"오류 발생: {e}")