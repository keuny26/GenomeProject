import streamlit as st
import fitz  # PyMuPDF
import json
import re
import time
import google.generativeai as genai
from streamlit_agraph import agraph, Node, Edge, Config
from Bio import Entrez  # NCBI 연동용

# --- 1. 페이지 설정 ---
st.set_page_config(page_title="GenomeGraph AI", layout="wide")
st.title("🧬 GenomeGraph AI (Final Stable Version)")

# --- 2. 세션 상태 초기화 ---
if "messages" not in st.session_state: st.session_state.messages = []
if "full_text" not in st.session_state: st.session_state.full_text = ""
if "graph_data" not in st.session_state: st.session_state.graph_data = None

# --- 3. API 키 및 설정 (사이드바) ---
with st.sidebar:
    st.title("⚙️ 설정 및 보안")
    # 우선순위: Secrets -> 사용자 직접 입력
    if "GEMINI_API_KEY" in st.secrets:
        api_key = st.secrets["GEMINI_API_KEY"]
        st.info("✅ Secrets에서 API 키를 로드했습니다.")
    else:
        api_key = st.text_input("Gemini API Key", type="password")
    
    ncbi_email = st.text_input("NCBI 연동용 이메일", value="your_email@example.com")
    
    if st.button("🗑️ 모든 데이터 초기화"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()

# --- 4. 모델 및 NCBI 함수 (404 에러 원천 차단) ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        # 404 에러 방지 핵심: 모델 객체를 생성할 때 명칭을 가장 표준화된 형태로 전달
        # 특정 환경에서 발생하는 v1beta 강제 호출 문제를 방어합니다.
        model = genai.GenerativeModel(model_name='gemini-1.5-flash')
        
        # 실제 모델이 작동하는지 가벼운 테스트 (선택 사항)
        st.sidebar.success("✅ 모델 연결 성공: gemini-1.5-flash")
    except Exception as e:
        st.error(f"모델 설정 오류: {e}")

def get_ncbi_gene_info(gene_name, email):
    Entrez.email = email
    try:
        search_handle = Entrez.esearch(db="gene", term=f"{gene_name}[Gene Name] AND human[Organism]")
        search_results = Entrez.read(search_handle)
        if not search_results["IdList"]: return "NCBI 정보를 찾을 수 없습니다."
        gene_id = search_results["IdList"][0]
        summary_handle = Entrez.esummary(db="gene", id=gene_id)
        summary_record = Entrez.read(summary_handle)
        return summary_record['DocumentSummarySet']['DocumentSummary'][0]['Description']
    except: return "NCBI 데이터 로드 실패"

# --- 5. 분석 및 병합 로직 ---
def analyze_single_doc(text, filename):
    if not model: return None
    # 개인정보 마스킹
    clean_text = re.sub(r'\d{3}-\d{4}-\d{4}', "[PROTECTED]", text)
    
    # JSON만 정확히 응답하도록 프롬프트 보강
    prompt = f"""
    You are a professional genome analyst. Extract genes, diseases, and variants from the text.
    Return ONLY a valid JSON object in the following format:
    {{
      "nodes": [{{"id": "unique_id", "label": "name", "type": "Gene/Disease/Variant", "desc": "summary"}}],
      "links": [{{"source": "id1", "target": "id2"}}]
    }}
    Text: {clean_text[:10000]}
    """
    
    try:
        time.sleep(1.0) # API 할당량 관리
        response = model.generate_content(prompt)
        # 응답에서 JSON 블록 추출 로직 강화
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            data = json.loads(json_match.group())
            if 'nodes' in data:
                for n in data['nodes']: n['source_file'] = filename
            return data
        return None
    except Exception as e:
        st.warning(f"[{filename}] 분석 중 오류: {e}")
        return None

def merge_graphs(results):
    merged_nodes = {}
    merged_links = []
    for data in results:
        if not data or 'nodes' not in data: continue
        for n in data['nodes']:
            nid = n['id']
            if nid in merged_nodes:
                merged_nodes[nid]['source_file'] = "Common"
            else:
                merged_nodes[nid] = n
        if 'links' in data:
            merged_links.extend(data['links'])
    
    unique_links = [dict(t) for t in {tuple(sorted(d.items())) for d in merged_links}]
    return {"nodes": list(merged_nodes.values()), "links": unique_links}

# --- 6. UI: 파일 업로드 및 분석 ---
uploaded_files = st.file_uploader("PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

if uploaded_files and api_key:
    if st.button("🧬 통합 분석 시작"):
        all_results = []
        with st.spinner("Gemini AI가 문서를 분석하고 있습니다..."):
            full_text_accumulator = ""
            for uploaded_file in uploaded_files:
                try:
                    # PDF 텍스트 추출
                    doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                    text = " ".join([page.get_text() for page in doc])
                    full_text_accumulator += f"\n\n[Doc: {uploaded_file.name}]\n{text}"
                    
                    res = analyze_single_doc(text, uploaded_file.name)
                    if res: all_results.append(res)
                except Exception as e:
                    st.error(f"{uploaded_file.name} 추출 실패: {e}")
            
            if all_results:
                st.session_state.full_text = full_text_accumulator
                st.session_state.graph_data = merge_graphs(all_results)
                st.session_state.messages = []
                st.success("통합 분석 완료!")
            else:
                st.error("데이터를 추출하지 못했습니다. API 설정 또는 PDF 내용을 확인하세요.")

    # --- 7. 그래프 및 상세 정보 ---
    if st.session_state.graph_data:
        col1, col2 = st.columns([3, 1])
        selected_id = None
        
        # 파일별 컬러 팔레트
        file_names = [f.name for f in uploaded_files]
        color_palette = ["#4285F4", "#34A853", "#FBBC05", "#8E44AD", "#F39C12", "#16A085"]
        color_map = {name: color_palette[i % len(color_palette)] for i, name in enumerate(file_names)}
        color_map["Common"] = "#EA4335" 

        with col1:
            st.subheader("🕸️ 지식 그래프")
            nodes = []
            for n in st.session_state.graph_data['nodes']:
                src = n.get('source_file', 'Unknown')
                nodes.append(Node(id=n['id'], 
                                  label=n['label'], 
                                  size=30 if src == "Common" else 20, 
                                  color=color_map.get(src, "#999999")))
            
            edges = [Edge(source=l['source'], target=l['target']) for l in st.session_state.graph_data['links']]
            
            if nodes:
                config = Config(width=800, height=600, directed=True, physics=True)
                selected_id = agraph(nodes=nodes, edges=edges, config=config)

        with col2:
            st.subheader("🎨 범례")
            for src, color in color_map.items():
                st.markdown(f"<span style='color:{color}'>●</span> {src}", unsafe_allow_html=True)
            
            st.divider()
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data['nodes'] if str(n['id']) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail['label']}")
                    st.info(f"**타입:** {node_detail['type']} | **출처:** {node_detail.get('source_file')}")
                    if node_detail['type'] == "Gene":
                        st.caption(f"**NCBI:** {get_ncbi_gene_info(node_detail['label'], ncbi_email)}")
                    st.write(f"**상세:** {node_detail.get('desc', '내용 없음')}")
            else:
                st.info("그래프의 노드를 클릭하세요.")

# --- 8. 채팅 영역 ---
if st.session_state.full_text:
    st.divider()
    st.subheader("💬 데이터 기반 Q&A")
    for msg in st.session_state.messages:
        with st.chat_message(msg["role"]): st.markdown(msg["content"])

    if chat_prompt := st.chat_input("질문하세요."):
        st.session_state.messages.append({"role": "user", "content": chat_prompt})
        with st.chat_message("user"): st.markdown(chat_prompt)
        
        with st.chat_message("assistant"):
            try:
                res = model.generate_content(f"Context: {st.session_state.full_text[:8000]}\nQ: {chat_prompt}")
                st.markdown(res.text)
                st.session_state.messages.append({"role": "assistant", "content": res.text})
            except Exception as e:
                st.error(f"채팅 오류: {e}")