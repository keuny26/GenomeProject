# Project/routers/document.py (Genomics Document Management & AI Extraction - 수정 완료)

from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from sqlalchemy.orm import Session
from typing import List, Optional, Dict, Any
import json
import os
from pydantic import BaseModel, Field

# ✅ AI 라이브러리 임포트
from google import genai
from google.genai import types

# ✅ 프로젝트 내부 모듈 임포트
import schemas
import crud 
from database import get_db

# ----------------------------------------------------------------------
# 1. LLM 응답을 위한 Pydantic 스키마 정의 (유지)
# ----------------------------------------------------------------------

class ExtractedGene(BaseModel):
    """추출된 유전자 정보"""
    symbol: str = Field(..., description="Gene symbol (e.g., CFTR, BRCA1)")
    name: str = Field(..., description="Full gene name.")
    location: str = Field(..., description="Chromosome location (e.g., 7q31.2).")

class ExtractedVariant(BaseModel):
    """추출된 변이 정보"""
    gene_symbol: str = Field(..., description="The symbol of the gene this variant belongs to.")
    name: str = Field(..., description="Variant name/ID (e.g., F508del, rs1801133).")
    variant_type: str = Field(..., description="Type of variant (e.g., SNP, Deletion, Missense).")

class ExtractedAssociation(BaseModel):
    """추출된 변이-질병 연관성 정보"""
    variant_name: str = Field(..., description="The Variant name/ID.")
    disease_name: str = Field(..., description="The Disease name.")
    p_value: Optional[float] = Field(None, description="Reported p-value for the association.")

class AI_ExtractionResult(BaseModel):
    """LLM이 최종적으로 반환해야 할 JSON 구조"""
    genes: List[ExtractedGene] = Field(..., description="Extracted Gene information list.")
    variants: List[ExtractedVariant] = Field(..., description="Extracted Variant information list.")
    associations: List[ExtractedAssociation] = Field(..., description="Extracted Association list.")
    
# ----------------------------------------------------------------------
# 2. AI 클라이언트 및 라우터 초기화 (유지)
# ----------------------------------------------------------------------

ai_client = None
try:
    # 환경 변수에서 API 키를 자동으로 로드하여 클라이언트 초기화
    ai_client = genai.Client() 
    print("✅ GenAI Client initialized.")
except Exception as e:
    # 실제 운영 환경에서는 키가 없으면 서버를 중단해야 하지만, 개발 편의상 print만 합니다.
    print(f"❌ GenAI Client initialization failed. Check your GEMINI_API_KEY in .env: {e}")

# ✅ APIRouter 객체 선언 (모든 @router.xxx 보다 앞에 와야 NameError 방지)
router = APIRouter(
    prefix="/documents",
    tags=["Documents (문서 관리 & AI 추출)"],
    responses={404: {"description": "Not found"}},
)

# ======================================================================
# 3. 문서 CRUD 엔드포인트 (수정: 목록 조회 스키마 변경)
# ======================================================================

# 🚨 수정됨: response_model을 schemas.DocumentListItem으로 변경합니다.
# 이는 schemas.py에서 정의한 간소화된 스키마로, 직렬화 오류를 방지합니다.
@router.get("/", response_model=List[schemas.DocumentListItem])
def read_documents(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """모든 문서 목록을 조회합니다."""
    documents = crud.get_documents(db, skip=skip, limit=limit)
    return documents

@router.get("/{document_id}", response_model=schemas.Document)
def read_document(document_id: int, db: Session = Depends(get_db)):
    """특정 ID의 문서를 조회합니다."""
    db_document = crud.get_document(db, document_id=document_id)
    if db_document is None:
        raise HTTPException(status_code=404, detail="Document not found")
    return db_document

# ======================================================================
# 4. 문서 업로드 및 AI 처리 (핵심 AI 기능 - 유지)
# ======================================================================

@router.post("/upload/", status_code=status.HTTP_201_CREATED)
async def upload_document_and_process(
    file: UploadFile = File(..., description="업로드할 임상 보고서 또는 논문 파일"),
    db: Session = Depends(get_db)
):
    """
    문서를 업로드하고, AI를 호출하여 핵심 엔티티와 관계를 추출한 후 DB에 저장합니다.
    (하나의 트랜잭션으로 문서 저장, AI 분석, 엔티티 저장을 처리합니다.)
    """
    if not ai_client:
        raise HTTPException(
            status_code=503, 
            detail="AI Client is not initialized. Check your GEMINI_API_KEY in .env file."
        )

    # 4.1. 파일 내용 읽기 및 문서 저장
    try:
        file_content_bytes = await file.read()
        document_text = file_content_bytes.decode('utf-8')
        file_name = file.filename
        
        # 4.1.1. 문서 자체를 DB에 저장 (트랜잭션 시작)
        db_document = crud.create_document(db=db, file_name=file_name, content=document_text)
        document_id = db_document.id
        
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=400, detail=f"파일 읽기/초기 저장 오류: {e}")

    # 4.2. LLM 시스템 명령어 및 설정 정의
    system_instruction = (
        "당신은 유전체 데이터 분석 전문가입니다. 제공된 문서에서 유전자(Gene), 변이(Variant), "
        "그리고 이들의 연관성(Association)에 대한 정보를 최대한 정확하게 추출해야 합니다. "
        "추출된 정보는 반드시 아래 JSON 스키마 형식에 맞춰서만 반환해야 합니다. "
        "추출할 정보가 없으면 해당 배열은 빈 리스트([])로 반환해야 합니다."
    )

    # 4.3. LLM 호출 (Structured Output 강제)
    try:
        response = ai_client.models.generate_content(
            model='gemini-2.5-flash',
            contents=[document_text],
            config=types.GenerateContentConfig(
                system_instruction=system_instruction,
                # LLM에게 JSON 형식으로 응답할 것을 명시
                response_mime_type="application/json",
                # Pydantic 스키마를 사용하여 JSON 구조를 강제
                response_schema=AI_ExtractionResult,
            )
        )
    except Exception as e:
        # LLM 호출 오류 시, 저장했던 문서까지 롤백
        print(f"LLM API Call Error: {e}")
        db.rollback()
        raise HTTPException(
            status_code=500, 
            detail=f"LLM API 호출 중 오류가 발생했습니다. (키 확인 필요): {e}"
        )

    # 4.4. LLM 응답 파싱 및 DB 저장 (트랜잭션 계속)
    try:
        raw_data = json.loads(response.text)
        
        saved_genes = {}
        # 4.4.1. 유전자(Gene) 저장 (get_or_create 사용)
        for gene_data in raw_data.get("genes", []):
            db_gene = crud.get_or_create_gene(db=db, **gene_data) 
            saved_genes[db_gene.symbol] = db_gene.id
            
        saved_variants = {}
        # 4.4.2. 변이(Variant) 저장 (문서 ID 연결)
        for variant_data in raw_data.get("variants", []):
            gene_symbol = variant_data.pop("gene_symbol")
            gene_id = saved_genes.get(gene_symbol) 
            
            if gene_id:
                db_variant = crud.create_variant(
                    db=db, 
                    variant=schemas.VariantCreate(
                        gene_id=gene_id, 
                        document_id=document_id, # 문서 ID 연결
                        **variant_data
                    )
                )
                saved_variants[db_variant.name] = db_variant.id

        # 4.4.3. 연관성(Association) 저장
        for assoc_data in raw_data.get("associations", []):
            variant_name = assoc_data.pop("variant_name")
            disease_name = assoc_data.pop("disease_name")
            
            db_disease = crud.get_or_create_disease(db=db, name=disease_name) 
            variant_id = saved_variants.get(variant_name)
            
            if variant_id and db_disease:
                crud.create_association(
                    db=db, 
                    association=schemas.AssociationCreate(
                        variant_id=variant_id, 
                        disease_id=db_disease.id, 
                        document_id=document_id, # 문서 ID 연결
                        p_value=assoc_data.get("p_value")
                    )
                )

        # 4.5. 최종 커밋 및 성공 응답
        db.commit() 
        return {
            "message": f"문서({file_name})가 성공적으로 처리되었고 AI 추출 결과가 데이터베이스에 저장되었습니다.",
            "extracted_counts": {
                "genes": len(raw_data.get("genes", [])),
                "variants": len(raw_data.get("variants", [])),
                "associations": len(raw_data.get("associations", [])),
            }
        }

    except Exception as e:
        db.rollback() 
        print(f"Database/Parsing Error: {e}")
        # 오류 발생 시 저장된 문서 자체도 롤백
        raise HTTPException(
            status_code=500, 
            detail=f"데이터베이스 저장/파싱 중 오류 발생 (롤백됨): {e}"
        )