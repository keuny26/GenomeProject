# Project/schemas.py (Genomics Schemas - 최종 통합 및 목록 최적화 버전)

from pydantic import BaseModel, Field, ConfigDict, model_validator
from typing import List, Optional, ForwardRef

# Pydantic v2에서 순환 참조를 위해 ForwardRef를 사용합니다.
Variant = ForwardRef("Variant")
Disease = ForwardRef("Disease")
Association = ForwardRef("Association")
Document = ForwardRef("Document") 

# ----------------------------------------------------------------------
# 0. Document (문서) 스키마
# ----------------------------------------------------------------------
class DocumentBase(BaseModel):
    name: str = Field(..., max_length=255, description="문서 제목 또는 파일 이름")
    content: str = Field(..., description="문서의 전체 텍스트 내용")

class DocumentCreate(DocumentBase):
    """문서 생성 요청 스키마 (업로드 시 사용)"""
    pass

# ✅ 새로 추가: 목록 조회용 스키마 (GET /documents/)
class DocumentListItem(BaseModel):
    """문서 목록 테이블에 표시할 간소화된 스키마"""
    id: int
    name: str
    
    # 프론트엔드 DocumentList.tsx에서 요구하는 필드 추가.
    # 모델에는 없지만, 직렬화 시 Pydantic이 기본값으로 채워 오류를 방지합니다.
    is_analyzed: bool = Field(True, description="문서 분석 완료 여부 (프론트엔드 호환성용)")

    model_config = ConfigDict(from_attributes=True)

class Document(DocumentBase):
    """문서 응답 스키마 (특정 ID 조회 시 사용)"""
    id: int
    
    # 💡 관계 설정: 이 문서에서 추출된 모든 연관성 목록을 포함합니다. (특정 ID 조회 시에만)
    associations: List[Association] = Field([], description="이 문서에서 추출된 연관성 목록")

    model_config = ConfigDict(from_attributes=True)
    
# ----------------------------------------------------------------------
# 1. Gene (유전자) 스키마
# ----------------------------------------------------------------------
class GeneBase(BaseModel):
    """유전자 모델의 기본 필드"""
    symbol: str = Field(..., max_length=50, description="유전자 심볼 (예: CFTR)")
    name: str = Field(..., max_length=255, description="유전자 전체 이름")
    location: Optional[str] = Field(None, max_length=50, description="염색체 위치 (예: 7q31.2)")

class GeneCreate(GeneBase):
    """유전자 생성 요청 스키마 (AI 추출 결과 저장)"""
    pass

class Gene(GeneBase):
    """유전자 응답 스키마 (ID 및 관계 포함)"""
    id: int
    variants: List[Variant] = Field([], description="이 유전자에 속한 변이 목록") 
    
    model_config = ConfigDict(from_attributes=True)
    
# ----------------------------------------------------------------------
# 2. Variant (변이) 스키마
# ----------------------------------------------------------------------

class VariantBase(BaseModel):
    """변이 모델의 기본 필드"""
    name: str = Field(..., max_length=50, description="변이 명칭 (예: F508del, rs1801133)")
    variant_type: Optional[str] = Field(None, max_length=50, description="변이 유형 (예: Deletion, Missense, SNP)")

class VariantCreate(VariantBase):
    """변이 생성 요청 스키마 (AI 추출 결과 저장)"""
    gene_id: int = Field(..., description="변이가 속한 유전자의 ID (FK)") 
    document_id: Optional[int] = Field(None, description="변이가 추출된 문서 ID") 

class Variant(VariantBase):
    """변이 응답 스키마"""
    id: int
    gene_id: int
    document_id: Optional[int] = None
    
    associations: List[Association] = Field([], description="이 변이와 관련된 연관성 목록")
    
    model_config = ConfigDict(from_attributes=True)

# ----------------------------------------------------------------------
# 3. Disease (질병) 스키마
# ----------------------------------------------------------------------

class DiseaseBase(BaseModel):
    """질병 모델의 기본 필드"""
    name: str = Field(..., max_length=255, description="질병 이름")
    description: Optional[str] = Field(None, description="질병에 대한 간단한 설명")

class DiseaseCreate(DiseaseBase):
    """질병 생성 요청 스키마 (AI 추출 결과 저장)"""
    pass
    
class Disease(DiseaseBase):
    """질병 응답 스키마"""
    id: int
    associations: List[Association] = Field([], description="이 질병과 관련된 연관성 목록")

    model_config = ConfigDict(from_attributes=True)

# ----------------------------------------------------------------------
# 4. Association (연관성) 스키마
# ----------------------------------------------------------------------

class AssociationBase(BaseModel):
    """연관성 모델의 기본 필드"""
    p_value: Optional[float] = Field(None, description="통계적 유의성 P-value")
    odds_ratio: Optional[float] = Field(None, description="오즈비 (Odds Ratio)")
    reference: Optional[str] = Field(None, max_length=255, description="참고 문헌 또는 출처")

class AssociationCreate(AssociationBase):
    """연관성 생성 요청 스키마 (AI 추출 결과 저장)"""
    variant_id: int = Field(..., description="관련된 변이 ID (FK)")
    disease_id: int = Field(..., description="관련된 질병 ID (FK)")
    document_id: Optional[int] = Field(None, description="연관성이 추출된 문서 ID")

class Association(AssociationBase):
    """연관성 응답 스키마"""
    id: int
    variant_id: int
    disease_id: int
    document_id: Optional[int] = None
    
    # 💡 관계 설정: Association 조회 시 연결된 Entity 정보를 포함합니다.
    variant: Optional[Variant] = Field(None, description="관련된 변이 객체") 
    disease: Optional[Disease] = Field(None, description="관련된 질병 객체")
    document: Optional[Document] = Field(None, description="추출된 출처 문서 객체")
    
    model_config = ConfigDict(from_attributes=True)


# ======================================================================
# 5. 모델 간의 순환 참조 해결 (Pydantic ForwardRef 처리)
# ======================================================================
# 모든 스키마 정의가 끝난 후 호출하여 순환 참조 문제를 해결합니다.
Document.model_rebuild()
Gene.model_rebuild()
Variant.model_rebuild()
Disease.model_rebuild()
Association.model_rebuild()