import streamlit as st
import fitz  # PyMuPDF
import json
import re
import google.generativeai as genai
from streamlit_agraph import agraph, Node, Edge, Config

# --- 페이지 설정 ---
st.set_page_config(page_title="GenomeGraph AI", layout="wide")
st.title("🧬 GenomeGraph AI (Streamlit)")

# --- API 키 설정 (보안 강화) ---
if "GEMINI_API_KEY" in st.secrets:
    api_key = st.secrets["GEMINI_API_KEY"]
else:
    st.sidebar.title("설정")
    api_key = st.sidebar.text_input("Gemini API Key를 입력하세요", type="password")

# --- 모델 초기화 (404 에러 방지 로직) ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        # 중요: 모델명 앞에 'models/'를 붙여 경로를 명확히 지정합니다.
        model = genai.GenerativeModel(model_name='models/gemini-1.5-flash')
    except Exception as e:
        st.error(f"API 설정 중 오류 발생: {e}")
else:
    st.warning("API 키가 설정되지 않았습니다. Secrets에 키를 추가하거나 사이드바에 직접 입력해주세요.")

# --- 분석 함수 ---
def analyze_text_with_ai(text):
    if not model:
        return None
    
    # f-string 내 중괄호 이스케이프({{, }}) 적용
    prompt = f"""
    당신은 유전체 데이터 분석가입니다. 아래 텍스트에서 유전자와 질환을 추출하여 JSON 그래프 데이터로 만드세요.
    결과는 반드시 다른 설명 없이 순수한 JSON 형식으로만 응답하세요.
    {{
      "nodes": [{{ "id": "ID", "label": "이름", "type": "gene/disease", "color": "#HEX", "desc": "설명" }}],
      "links": [{{ "source": "ID", "target": "ID" }}]
    }}
    텍스트: {text[:8000]}
    """
    
    try:
        response = model.generate_content(prompt)
        res_text = response.text
        
        # JSON 추출 로직 (Markdown 및 공백 제거)
        clean_json = re.sub(r'```json|```', '', res_text).strip()
        return json.loads(clean_json)
    except Exception as e:
        st.error(f"AI 응답 처리 중 오류: {e}")
        return None

# --- 메인 UI ---
uploaded_file = st.file_uploader("PDF 보고서를 업로드하세요", type="pdf")

if uploaded_file and api_key:
    with st.spinner("AI가 유전체 데이터를 분석 중입니다..."):
        try:
            # 1. PDF 읽기
            content = uploaded_file.read()
            doc = fitz.open(stream=content, filetype="pdf")
            full_text = "".join([page.get_text() for page in doc])
            
            # 2. AI 분석
            graph_data = analyze_text_with_ai(full_text)
            
            if graph_data:
                # 3. 데이터 시각화 요소 생성
                nodes = []
                for n in graph_data.get('nodes', []):
                    # 타입에 따른 색상 구분
                    default_color = '#4285F4' if n.get('type') == 'gene' else '#EA4335'
                    nodes.append(Node(id=n['id'], label=n['label'], size=20, color=n.get('color', default_color)))
                
                edges = []
                for l in graph_data.get('links', []):
                    edges.append(Edge(source=l['source'], target=l['target']))

                # 4. 출력 레이아웃
                st.subheader("🧬 분석 결과 지식 그래프")
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    config = Config(width=900, height=600, directed=True, physics=True)
                    selected_id = agraph(nodes=nodes, edges=edges, config=config)
                
                with col2:
                    st.markdown("### 🔍 상세 정보")
                    if selected_id:
                        node_detail = next((n for n in graph_data['nodes'] if n['id'] == selected_id), None)
                        if node_detail:
                            st.success(f"**명칭:** {node_detail['label']}")
                            st.info(f"**설명:** {node_detail.get('desc', '설명 없음')}")
                    else:
                        st.write("노드를 클릭하면 상세 내용이 표시됩니다.")
        except Exception as e:
            st.error(f"분석 중 오류 발생: {e}")