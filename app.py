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
st.title("🧬 GenomeGraph AI (Universal Compatibility)")

# --- 2. 세션 상태 초기화 ---
if "messages" not in st.session_state: st.session_state.messages = []
if "full_text" not in st.session_state: st.session_state.full_text = ""
if "graph_data" not in st.session_state: st.session_state.graph_data = None
if "active_model_name" not in st.session_state: st.session_state.active_model_name = None

# --- 3. API 키 및 설정 (사이드바) ---
with st.sidebar:
    st.title("⚙️ 설정 및 보안")
    if "GEMINI_API_KEY" in st.secrets:
        api_key = st.secrets["GEMINI_API_KEY"]
        st.info("✅ Secrets에서 API 키를 로드했습니다.")
    else:
        api_key = st.text_input("Gemini API Key", type="password")
    
    ncbi_email = st.text_input("NCBI 연동용 이메일", value="your_email@example.com")
    
    if st.button("🗑️ 모든 데이터 초기화"):
        for key in list(st.session_state.keys()): del st.session_state[key]
        st.rerun()

# --- 4. 모델 설정 (진단 및 강제 우회 통합) ---
model = None
if api_key:
    try:
        genai.configure(api_key=api_key)
        
        # [진단] 현재 이 API 키로 접근 가능한 모든 모델 리스트 출력
        # 404 에러의 원인을 찾기 위해 전체 목록을 먼저 스캔합니다.
        all_models = [m.name for m in genai.list_models()]
        available_models = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]
        
        if not available_models:
            st.sidebar.error("❌ 접근 가능한 모델이 0개입니다.")
            st.sidebar.info(f"전체 모델 리스트: {all_models}")
            
            # [최후의 수단] 리스트에 없어도 강제로 이름을 지정하여 시도 (SDK 버그 우회)
            st.sidebar.warning("⚠️ 강제 연결 모드 진입: gemini-1.5-flash")
            model = genai.GenerativeModel(model_name='gemini-1.5-flash')
            st.session_state.active_model_name = "gemini-1.5-flash (Forced)"
        else:
            # 정상적인 경우 우선순위에 따라 선택 (1.5 Flash 선호)
            target = next((m for m in available_models if "gemini-1.5-flash" in m), available_models[0])
            model = genai.GenerativeModel(model_name=target)
            st.session_state.active_model_name = target
            st.sidebar.success(f"✅ 연결 성공: {target}")
            
    except Exception as e:
        st.error(f"❌ API 접근 오류 발생: {e}")
        if "API_KEY_INVALID" in str(e):
            st.info("💡 API 키가 유효하지 않습니다. AI Studio에서 키를 새로 생성해 보세요.")
        elif "PERMISSION_DENIED" in str(e):
            st.info("💡 권한 거부: Google Cloud Console에서 'Generative Language API'가 활성화되었는지 확인하세요.")

# --- 5. 분석 및 병합 로직 (기능 유지) ---
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

def analyze_single_doc(text, filename):
    if not model: return None
    clean_text = re.sub(r'\d{3}-\d{4}-\d{4}', "[PROTECTED]", text)
    prompt = f"""
    당신은 전문 유전체 분석가입니다. 다음 텍스트에서 유전자(Gene), 질환/증상(Disease), 변이(Variant) 관계를 추출하여 JSON으로만 응답하세요.
    - 출력 형식: {{"nodes": [{{"id": "ID", "label": "이름", "type": "Gene/Disease/Variant", "desc": "설명"}}], "links": [{{"source": "ID", "target": "ID"}}]}}
    텍스트: {clean_text[:10000]}
    """
    try:
        time.sleep(1.2)
        response = model.generate_content(prompt)
        json_match = re.search(r'\{.*\}', response.text, re.DOTALL)
        if json_match:
            data = json.loads(json_match.group())
            if 'nodes' in data:
                for n in data['nodes']: n['source_file'] = filename
            return data
    except Exception as e:
        st.warning(f"[{filename}] 분석 실패: {e}")
    return None

def merge_graphs(results):
    merged_nodes, merged_links = {}, []
    for data in results:
        if not data or 'nodes' not in data: continue
        for n in data['nodes']:
            nid = n['id']
            if nid in merged_nodes: merged_nodes[nid]['source_file'] = "Common"
            else: merged_nodes[nid] = n
        if 'links' in data: merged_links.extend(data['links'])
    unique_links = [dict(t) for t in {tuple(sorted(d.items())) for d in merged_links}]
    return {"nodes": list(merged_nodes.values()), "links": unique_links}

# --- 6. UI: 파일 업로드 및 분석 ---
uploaded_files = st.file_uploader("PDF 보고서 업로드 (다중 선택 가능)", type="pdf", accept_multiple_files=True)

if uploaded_files and api_key:
    if st.button("🧬 통합 분석 시작"):
        all_results = []
        with st.spinner(f"{st.session_state.active_model_name} 엔진으로 분석 중..."):
            full_txt = ""
            for f in uploaded_files:
                doc = fitz.open(stream=f.read(), filetype="pdf")
                text = " ".join([page.get_text() for page in doc])
                full_txt += f"\n\n[Doc: {f.name}]\n{text}"
                res = analyze_single_doc(text, f.name)
                if res: all_results.append(res)
            
            if all_results:
                st.session_state.full_text = full_txt
                st.session_state.graph_data = merge_graphs(all_results)
                st.success("분석 완료!")

    # --- 7. 그래프 필터링 및 시직화 ---
    if st.session_state.graph_data:
        st.sidebar.divider()
        st.sidebar.subheader("🔍 그래프 필터 및 검색")
        all_types = list(set([n.get('type', 'Unknown') for n in st.session_state.graph_data['nodes']]))
        selected_types = st.sidebar.multiselect("표시할 타입", all_types, default=all_types)
        search_query = st.sidebar.text_input("🎯 노드 검색 (이름)")

        col1, col2 = st.columns([3, 1])
        selected_id = None
        
        file_names = [f.name for f in uploaded_files]
        color_palette = ["#4285F4", "#34A853", "#FBBC05", "#8E44AD", "#F39C12", "#16A085"]
        color_map = {name: color_palette[i % len(color_palette)] for i, name in enumerate(file_names)}
        color_map["Common"] = "#EA4335" 

        with col1:
            f_nodes = []
            f_node_ids = set()
            for n in st.session_state.graph_data['nodes']:
                if n.get('type') in selected_types and search_query.lower() in n.get('label', '').lower():
                    src = n.get('source_file', 'Unknown')
                    f_nodes.append(Node(id=n['id'], label=n['label'], size=25 if src=="Common" else 20, color=color_map.get(src, "#999999")))
                    f_node_ids.add(n['id'])
            f_edges = [Edge(source=l['source'], target=l['target']) for l in st.session_state.graph_data['links'] 
                       if l['source'] in f_node_ids and l['target'] in f_node_ids]

            if f_nodes:
                selected_id = agraph(nodes=f_nodes, edges=f_edges, config=Config(width=900, height=600, directed=True, physics=True))

        with col2:
            st.markdown("### 🎨 범례 및 상세")
            for src, color in color_map.items():
                st.markdown(f"<span style='color:{color}'>●</span> {src}", unsafe_allow_html=True)
            st.divider()
            
            if selected_id:
                node = next((n for n in st.session_state.graph_data['nodes'] if str(n['id']) == str(selected_id)), None)
                if node:
                    st.success(f"**명칭:** {node['label']}")
                    st.info(f"**타입:** {node['type']} | **출처:** {node.get('source_file')}")
                    if node['type'] == "Gene":
                        with st.spinner("NCBI 확인 중..."):
                            st.caption(f"**NCBI:** {get_ncbi_gene_info(node['label'], ncbi_email)}")
                        st.link_button("🧬 NCBI 상세보기", f"https://www.ncbi.nlm.nih.gov/gene/?term={node['label']}")
                    st.write(f"**상세 설명:**\n{node.get('desc', '내용 없음')}")
            else:
                st.info("그래프의 노드를 선택하세요.")

# --- 8. 채팅 영역 ---
if st.session_state.full_text:
    st.divider()
    st.subheader("💬 데이터 기반 Q&A")
    for msg in st.session_state.messages:
        with st.chat_message(msg["role"]): st.markdown(msg["content"])
    if chat_prompt := st.chat_input("질문하세요."):
        st.session_state.messages.append({"role": "user", "content": chat_prompt})
        with st.chat_message("user"): st.markdown(chat_prompt)
        try:
            res = model.generate_content(f"Context: {st.session_state.full_text[:8000]}\nQ: {chat_prompt}")
            with st.chat_message("assistant"):
                st.markdown(res.text)
                st.session_state.messages.append({"role": "assistant", "content": res.text})
        except Exception as e: st.error(f"채팅 오류: {e}")