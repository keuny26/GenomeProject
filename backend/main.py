import fitz  # PyMuPDF
import json
import google.generativeai as genai
from fastapi import FastAPI, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

app = FastAPI()

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ✅ Gemini 설정 (발급받은 키를 여기에 넣으세요)
genai.configure(api_key="YOUR_GEMINI_API_KEY")
model = genai.GenerativeModel('gemini-1.5-flash')

class ChatQuery(BaseModel):
    query: str

async def analyze_text_with_ai(text: str):
    """AI에게 텍스트를 보내 구조화된 JSON 그래프 데이터를 받아옵니다."""
    prompt = f"""
    당신은 숙련된 유전체학 데이터 분석가입니다. 
    다음 텍스트에서 언급된 유전자(Gene)와 질환(Disease)을 추출하고 그들 사이의 연관성을 지식 그래프 형태로 요약하세요.
    
    결과는 반드시 아래의 JSON 형식으로만 응답하세요:
    {{
      "nodes": [
        {{"id": "unique_id", "label": "이름", "type": "gene 또는 disease", "color": "gene이면 #4285F4, disease면 #EA4335", "description": "상세 설명"}},
        ...
      ],
      "links": [
        {{"source": "node_id_1", "target": "node_id_2"}},
        ...
      ]
    }}

    텍스트 내용:
    {text[:8000]}  # 텍스트가 너무 길 경우를 대비해 슬라이싱
    """
    
    try:
        response = model.generate_content(prompt)
        # JSON 부분만 추출 (AI가 가끔 ```json ... ``` 코드를 붙이기 때문)
        json_str = response.text.replace('```json', '').replace('```', '').strip()
        return json.loads(json_str)
    except Exception as e:
        print(f"AI 분석 에ers: {e}")
        return None

@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    content = await file.read()
    try:
        # 1. PDF 텍스트 추출
        doc = fitz.open(stream=content, filetype="pdf")
        text = "".join([page.get_text() for page in doc])
        
        # 2. AI를 이용한 그래프 데이터 생성
        graph_data = await analyze_text_with_ai(text)
        
        if not graph_data:
            return {"graph": {"nodes": [], "links": []}}

        # 3. 문서 노드 추가 (중심점 역할을 위해 수동 추가)
        doc_node = {
            "id": f"doc_{file.filename}", 
            "label": file.filename, 
            "type": "doc", 
            "color": "#1a73e8", 
            "description": "업로드된 원본 문서"
        }
        graph_data["nodes"].append(doc_node)
        
        # 모든 첫 노드들을 문서 노드에 연결 (구조적 통일성)
        for node in graph_data["nodes"]:
            if node["type"] != "doc":
                graph_data["links"].append({"source": doc_node["id"], "target": node["id"]})

        return {"graph": graph_data}
    except Exception as e:
        print(f"Error: {e}")
        return {"graph": {"nodes": [], "links": []}}

@app.post("/chat")
async def chat_with_ai(data: ChatQuery):
    # 나중에 여기에 PDF 내용을 기억하는 RAG 기능을 추가할 수 있습니다.
    return {"reply": "AI 분석이 완료되었습니다. 그래프의 노드를 클릭하여 관계를 확인하세요."}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)