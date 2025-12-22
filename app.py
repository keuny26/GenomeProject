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
    if "GEMINI_API_KEY" in st.secrets:
        api_key = st.secrets["GEMINI_API_KEY"]
        st.info("✅ Secrets에서 API 키를 로드했습니다.")
    else:
        api_key = st.text_input("Gemini API Key", type="password", help="키는 세션 종료 시 삭제됩니다.")
    
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
        # 안정성을 위해 최신 flash 모델 지정
        model = genai.GenerativeModel('gemini-1.5-flash')
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
    # 보안: 개인정보 마스킹
    clean_text = re.sub(r'\d{3}-\d{4}-\d{4}', "[PROTECTED]", text)
    prompt = f"""
    당신은 전문 유전체 분석가입니다. 다음 텍스트에서 유전자, 질환(증상 포함), 변이 관계를 추출하여 JSON으로 응답하세요.
    - 증상/질병은 반드시 'type': 'Disease'로 분류하세요.
    - 출력 형식: {{"nodes": [{{"id": "ID", "label": "이름", "type": "Gene/Disease/Variant/Drug", "desc": "설명"}}], "links": [{{"source": "ID", "target": "ID"}}]}}
    텍스트: {clean_text[:10000]}
    """
    try:
        time.sleep(1) # API 할당량 조절
        response = model.generate_content(prompt)
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            data = json.loads(json_match.group())
            # 노드별로 출처 파일 정보 기록
            for n in data['nodes']: n['source_file'] = filename
            return data
        return None
    except Exception as e:
        st.error(f"{filename} 분석 오류: {e}")
        return None

def merge_graphs(results):
    merged_nodes = {}
    merged_links = []
    for data in results:
        if not data: continue
        for n in data['nodes']:
            nid = n['id']
            if nid in merged_nodes:
                # 여러 파일에서 발견된 노드는 'Common'으로 표시
                merged_nodes[nid]['source_file'] = "Common"
            else:
                merged_nodes[nid] = n
        merged_links.extend(data['links'])
    
    # 중복 링크 제거
    unique_links = [dict(t) for t in {tuple(sorted(d.items())) for d in merged_links}]
    return {"nodes": list(merged_nodes.values()), "links": unique_links}

# --- 6. UI: 파일 업로드 및 분석 ---
uploaded_files = st.file_uploader("PDF 보고서들을 업로드하세요", type="pdf", accept_multiple_files=True)

if uploaded_files and api_key:
    if st.button("🧬 파일 통합 분석 시작 (누락 방지 모드)"):
        all_results = []
        with st.spinner("문서별로 개별 분석을 진행 중입니다..."):
            full_text_accumulator = ""
            for uploaded_file in uploaded_files:
                doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
                text = " ".join([page.get_text() for page in doc])
                full_text_accumulator += f"\n\n[Document: {uploaded_file.name}]\n{text}"
                
                # 문서별 개별 분석 실행
                res = analyze_single_doc(text, uploaded_file.name)
                all_results.append(res)
            
            st.session_state.full_text = full_text_accumulator
            st.session_state.graph_data = merge_graphs(all_results)
            st.session_state.messages = []
            st.success("모든 문서 분석 및 통합 완료!")

    # --- 7. 그래프 시각화 및 필터 영역 ---
    if st.session_state.graph_data:
        st.sidebar.divider()
        st.sidebar.subheader("🔍 필터링")
        all_types = list(set([n.get('type', 'Unknown') for n in st.session_state.graph_data['nodes']]))
        selected_types = st.sidebar.multiselect("노드 타입", all_types, default=all_types)
        search_query = st.sidebar.text_input("🎯 노드 검색")

        col1, col2 = st.columns([3, 1])
        
        # NameError 방지를 위해 selected_id 초기화
        selected_id = None
        
        # 파일별 색상 설정
        file_names = [f.name for f in uploaded_files]
        color_palette = ["#4285F4", "#34A853", "#FBBC05", "#8E44AD", "#F39C12", "#16A085"]
        color_map = {name: color_palette[i % len(color_palette)] for i, name in enumerate(file_names)}
        color_map["Common"] = "#EA4335" 

        with col1:
            filtered_nodes = []
            filtered_node_ids = set()
            for n in st.session_state.graph_data['nodes']:
                type_match = n.get('type') in selected_types
                search_match = search_query.lower() in n.get('label', '').lower()
                
                if type_match and search_match:
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
                config = Config(width=900, height=600, directed=True, physics=True, fit_view=True, panAndZoom=True)
                # 여기서 클릭된 노드 ID를 가져옴
                selected_id = agraph(nodes=filtered_nodes, edges=filtered_edges, config=config)

        with col2:
            st.markdown("### 🎨 범례 및 상세 정보")
            for src, color in color_map.items():
                st.markdown(f"<span style='color:{color}'>●</span> **{src}**", unsafe_allow_html=True)
            
            st.divider()
            # selected_id가 None이 아닐 때만 상세 정보 출력 (NameError 방지)
            if selected_id:
                node_detail = next((n for n in st.session_state.graph_data['nodes'] if str(n['id']) == str(selected_id)), None)
                if node_detail:
                    st.success(f"**명칭:** {node_detail['label']} ({node_detail['type']})")
                    st.info(f"**출처:** {node_detail.get('source_file', 'N/A')}")
                    
                    if node_detail['type'] == "Gene":
                        with st.spinner("NCBI 확인 중..."):
                            ncbi_info = get_ncbi_gene_info(node_detail['label'], ncbi_email)
                            st.caption(f"**NCBI Summary:** {ncbi_info}")
                    
                    st.link_button("🧬 NCBI 바로가기", f"https://www.ncbi.nlm.nih.gov/gene/?term={node_detail['label']}")
                    st.write(f"**AI 분석 상세:**\n{node_detail.get('desc', '설명 없음')}")
            else:
                st.write("그래프의 노드를 클릭하면 상세 정보를 볼 수 있습니다.")

# --- 8. 채팅 영역 ---
st.divider()
st.subheader("💬 데이터 보안 채팅")
for message in st.session_state.messages:
    with st.chat_message(message["role"]): st.markdown(message["content"])

if chat_prompt := st.chat_input("데이터에 대해 질문하세요."):
    st.session_state.messages.append({"role": "user", "content": chat_prompt})
    with st.chat_message("user"): st.markdown(chat_prompt)
    with st.chat_message("assistant"):
        try:
            # 전체 텍스트 기반 답변 (8000자 제한)
            res = model.generate_content(f"Context: {st.session_state.full_text[:8000]}\nQuestion: {chat_prompt}")
            st.markdown(res.text)
            st.session_state.messages.append({"role": "assistant", "content": res.text})
        except Exception as e: st.error(f"채팅 오류: {e}")