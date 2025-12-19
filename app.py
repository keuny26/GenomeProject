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

# --- 모델 초기화 (404 에러 및 버전 문제 해결) ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        
        # 내 키로 사용 가능한 모델 목록을 확인하여 자동으로 매칭합니다.
        available_models = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]
        
        if 'models/gemini-1.5-flash' in available_models:
            model = genai.GenerativeModel('gemini-1.5-flash')
        elif 'models/gemini-1.5-pro' in available_models:
            model = genai.GenerativeModel('gemini-1.5-pro')
        elif available_models:
            # 리스트에 있는 첫 번째 가용 모델 선택
            model = genai.GenerativeModel(available_models[0].replace('models/', ''))
        
        if model:
            st.sidebar.success(f"연결됨: {model.model_name}")
        else:
            st.error("사용 가능한 Gemini 모델을 찾을 수 없습니다.")
            
    except Exception as e:
        st.error(f"API 설정 중 오류 발생: {e}")
else:
    st.warning("API 키가 설정되지 않았습니다. Secrets에 키를 추가하거나 사이드바에 직접 입력해주세요.")

# --- 분석 함수 ---
def analyze_text_with_ai(text):
    if not model:
        return None
    
    safe_text = text[:15000]
    prompt = f"""
    당신은 유전체 데이터 분석가입니다. 아래 텍스트에서 유전자와 질환의 관계를 추출하여 JSON 그래프 데이터로 만드세요.
    반드시 JSON 형식으로만 답변하세요.

    {{
      "nodes": [{{ "id": "ID", "label": "이름", "type": "gene/disease", "desc": "설명" }}],
      "links": [{{ "source": "ID", "target": "ID" }}]
    }}

    텍스트: {safe_text}
    """
    
    try:
        response = model.generate_content(prompt)
        res_text = response.text
        
        json_match = re.search(r'\{.*\}', res_text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group())
        return None
    except Exception as e:
        st.error(f"AI 분석 중 오류: {e}")
        return None

# --- 메인 UI ---
uploaded_file = st.file_uploader("PDF 보고서를 업로드하세요", type="pdf")

if uploaded_file and api_key:
    with st.spinner("AI가 유전체 데이터를 분석 중입니다..."):
        try:
            content = uploaded_file.read()
            doc = fitz.open(stream=content, filetype="pdf")
            full_text = " ".join([page.get_text() for page in doc])
            
            if not full_text.strip():
                st.error("PDF에서 텍스트를 추출할 수 없습니다.")
            else:
                graph_data = analyze_text_with_ai(full_text)
                
                if graph_data and 'nodes' in graph_data:
                    nodes = []
                    for n in graph_data.get('nodes', []):
                        color = '#4285F4' if n.get('type') == 'gene' else '#EA4335'
                        nodes.append(Node(id=n['id'], label=n['label'], size=25, color=color))
                    
                    edges = []
                    for l in graph_data.get('links', []):
                        edges.append(Edge(source=l['source'], target=l['target']))

                    st.subheader("🧬 분석 결과 지식 그래프")
                    col1, col2 = st.columns([3, 1])
                    
                    with col1:
                        config = Config(width=800, height=600, directed=True, physics=True)
                        selected_id = agraph(nodes=nodes, edges=edges, config=config)
                    
                    with col2:
                        st.markdown("### 🔍 상세 정보")
                        if selected_id:
                            node_detail = next((n for n in graph_data['nodes'] if n['id'] == selected_id), None)
                            if node_detail:
                                st.success(f"**명칭:** {node_detail['label']}")
                                st.info(f"**설명:** {node_detail.get('desc', '설명 없음')}")
                        else:
                            st.write("노드를 클릭하세요.")
        except Exception as e:
            st.error(f"프로세싱 중 오류: {e}")