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
    # GitHub Secrets 사용 권장
    if "GEMINI_API_KEY" in st.secrets:
        api_key = st.secrets["GEMINI_API_KEY"]
        st.info("✅ Secrets에서 API 키를 로드했습니다.")
    else:
        api_key = st.text_input("Gemini API Key", type="password", help="GitHub Secrets 사용을 권장합니다.")
    
    ncbi_email = st.text_input("NCBI 연동용 이메일", value="your_email@example.com")
    
    if st.button("🗑️ 모든 데이터 초기화"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()

# --- 4. 모델 및 NCBI 함수 ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        # 404 에러 방지: 모델 이름만 명확히 전달
        model = genai.GenerativeModel(model_name='gemini-1.5-flash')
        st.sidebar.success("모델 연결됨: gemini-1.5-flash")
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
    # 보안: 개인정보 마스킹 (전화번호 형태 등)
    clean_text = re.sub(r'\d{3}-\d{4}-\d{4}', "[PROTECTED]", text)
    
    prompt = f"""
    당신은 전문 유전체 분석가입니다. 다음 텍스트에서 유전자(Gene), 질환/증상(Disease), 변이(Variant) 관계를 추출하여 JSON으로만 응답하세요.
    - 증상/질병은 반드시 'type': 'Disease'로 분류하세요.
    - 출력 형식: {{"nodes": [{{"id": "ID", "label": "이름", "type": "Gene/Disease/Variant", "desc": "설명"}}], "links": [{{"source": "ID", "target": "ID"}}]}}
    텍스트: {clean_text[:12000]}
    """
    
    try:
        time.sleep(1) # API 할당량 관리
        response = model.generate_content(prompt)
        # JSON 블록만 추출
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            data = json.loads(json_match.group())
            # 노드별 출처 정보 기록
            for n in data.get('nodes', []): n['source_file'] = filename
            return data
        return None
    except Exception as e:
        st.error(f"{filename} 분석 중 오류 발생: {e}")
        return None

def merge_graphs(results):
    merged_nodes = {}
    merged_links = []
    for data in results:
        if not data: continue
        for n in data.get('nodes', []):
            nid = n['id']
            if nid in merged_nodes:
                merged_nodes[nid]['source_file'] = "Common"
            else:
                merged_nodes[nid] = n
        merged_links.extend(data.get('links', []))
    
    # 중복 링크 제거
    unique_links = [dict(t) for t in {tuple(sorted(d.items())) for d in merged_links}]
    return {"nodes": list(merged_nodes.values()), "links": unique_links}

# --- 6. UI: 파일 업로드 및 분석 ---
uploaded_files = st.file_uploader("분석할 PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

if uploaded_files and api_key:
    if st.button("🧬 통합 분석 시작 (Multi-Doc Mode)"):
        all_results = []
        with st.spinner("문서별 정밀 분석 진행 중..."):
            full_text_accumulator = ""
            for uploaded_file in uploaded_files:
                doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                text = " ".join([page.get_text() for page in doc])
                full_text_accumulator += f"\n\n[Document: {uploaded_file.name}]\n{text}"
                
                res = analyze_single_doc(text, uploaded_file.name)
                all_results.append(res)
            
            st.session_state.full_text = full_text_accumulator
            st.session_state.graph_data = merge_graphs(all_results)
            st.session_state.messages = []
            st.success("통합 분석 완료!")

    # --- 7. 그래프 시각화 영역 ---
    if st.session_state.graph_data:
        st.sidebar.divider()
        st.sidebar.subheader("🔍 필터링")
        all_types = list(set([n.get('type', 'Unknown') for n in st.session_state.graph_data['nodes']]))
        selected_types = st.sidebar.multiselect("표시할 타입", all_types, default=all_types)
        search_query = st.sidebar.text_input("🎯 노드 검색")

        col1, col2 = st.columns([3, 1])
        
        # NameError 방지: selected_id 초기화
        selected_id = None
        
        # 파일별 컬러 맵핑
        file_names = [f.name for f in uploaded_files]
        color_palette = ["#4285F4", "#34A853", "#FBBC05", "#8E44AD", "#F39C12", "#16A085"]
        color_map = {name: color_palette[i % len(color_palette)] for i, name in enumerate(file_names)}
        color_map["Common"] = "#EA4335" 

        with col1:
            filtered_nodes = []
            filtered_node_ids = set()
            for n in st.session_state.graph_data['nodes']:
                if n.get('type') in selected_types and search_query.lower() in n.get('label', '').lower():
                    src = n.get('source_file', 'Unknown')
                    n_color = color_map.get(src, "#999999")
                    is_common = src == "Common"
                    filtered_nodes.append(Node(id=n['id'], 
                                               label=f"⭐ {n['label']}" if is_common else n['label'], 
                                               size=35 if is_common else 25, 
                                               color=n_color))
                    filtered_node_ids.add(n['id'])
            
            filtered_edges = [Edge(source=l['source'], target=l['target']) 
                              for l in st.session_state.graph_data['links'] 
                              if l['source'] in filtered_node_ids and l['target'] in filtered_node_ids]

            if filtered_nodes:
                config = Config(width=900, height=600, directed=True, physics=True, fit_view=True)
                selected_id = agraph(nodes=filtered_nodes, edges=filtered_edges, config=config)

        with col2:
            st.markdown("### 🎨 범례")
            for src, color in color_map.items():
                st.markdown(f"<span style='color:{color}'>●</span> **{src}**", unsafe_allow_html=True)
            
            st.divider()
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data['nodes'] if str(n['id']) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**이름:** {node_detail['label']}")
                    st.info(f"**타입:** {node_detail['type']} | **출처:** {node_detail.get('source_file')}")
                    
                    if node_detail['type'] == "Gene":
                        with st.spinner("NCBI 데이터 검색 중..."):
                            ncbi_info = get_ncbi_gene_info(node_detail['label'], ncbi_email)
                            st.caption(f"**NCBI Summary:** {ncbi_info}")
                    
                    st.link_button("🧬 NCBI 상세보기", f"https://www.ncbi.nlm.nih.gov/gene/?term={node_detail['label']}")
                    st.write(f"**상세 설명:**\n{node_detail.get('desc', '내용 없음')}")
            else:
                st.info("그래프 노드를 클릭하면 상세 정보가 표시됩니다.")

# --- 8. 채팅 영역 ---
if st.session_state.full_text:
    st.divider()
    st.subheader("💬 분석 데이터 기반 Q&A")
    for message in st.session_state.messages:
        with st.chat_message(message["role"]): st.markdown(message["content"])

    if chat_prompt := st.chat_input("이 유전체 데이터들에 대해 질문하세요."):
        st.session_state.messages.append({"role": "user", "content": chat_prompt})
        with st.chat_message("user"): st.markdown(chat_prompt)
        
        with st.chat_message("assistant"):
            try:
                # 404 방지를 위해 초기화된 model 객체 사용
                res = model.generate_content(f"Context: {st.session_state.full_text[:10000]}\nQuestion: {chat_prompt}")
                st.markdown(res.text)
                st.session_state.messages.append({"role": "assistant", "content": res.text})
            except Exception as e:
                st.error(f"채팅 응답 오류: {e}")