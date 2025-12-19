# Project/routers/gene.py

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import List, Optional

# 로컬 모듈 임포트
import crud, schemas
from database import get_db
from models import GeneModel # 모델 임포트는 필요 없지만, 타입 힌트를 위해 포함 가능

# 라우터 인스턴스 생성
router = APIRouter(
    prefix="/genes",
    tags=["Genes (유전자)"],
)

# ======================================================================
# 1. 유전자 목록 조회 (Read All)
# ======================================================================
@router.get("/", response_model=List[schemas.Gene])
def read_genes(
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """
    모든 유전자 목록을 조회합니다.
    """
    genes = crud.get_genes(db, skip=skip, limit=limit)
    return genes

# ======================================================================
# 2. 유전자 생성 (Create)
# ======================================================================
@router.post("/", response_model=schemas.Gene, status_code=status.HTTP_201_CREATED)
def create_new_gene(gene: schemas.GeneCreate, db: Session = Depends(get_db)):
    """
    새로운 유전자를 생성합니다. (유전자 심볼은 중복될 수 없습니다.)
    """
    # 유전자 심볼 중복 체크
    db_gene = crud.get_gene_by_symbol(db, symbol=gene.symbol)
    if db_gene:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Gene with symbol '{gene.symbol}' already exists."
        )
    
    # crud 함수 호출 시 IntegrityError가 발생할 수도 있으므로 try-except 사용 (create_gene에 IntegrityError 처리 포함됨)
    try:
        return crud.create_gene(db=db, gene=gene)
    except ValueError as e:
        # crud.create_gene에서 발생시킨 ValueError 처리
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, 
            detail=str(e)
        )

# ======================================================================
# 3. 특정 유전자 조회 (Read One)
# ======================================================================
@router.get("/{gene_id}", response_model=schemas.Gene)
def read_gene(gene_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 유전자를 조회합니다."""
    # 관계형 데이터(variants)도 함께 로드하기 위해 .options(selectinload(GeneModel.variants))를 사용해야 하지만,
    # Pydantic 스키마 (schemas.Gene)에 variants 리스트가 정의되어 있으므로, 
    # 모델에 관계가 설정되어 있다면 자동으로 로드됩니다 (Lazy Loading 또는 모델 설정에 따라).
    db_gene = crud.get_gene(db, gene_id=gene_id)
    if db_gene is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"유전자 ID {gene_id}를 찾을 수 없습니다."
        )
    return db_gene

# ======================================================================
# 4. 유전자 수정 (Update)
# ======================================================================
@router.put("/{gene_id}", response_model=schemas.Gene)
def update_existing_gene(gene_id: int, gene_in: schemas.GeneCreate, db: Session = Depends(get_db)):
    """특정 ID를 가진 유전자 정보를 수정합니다."""
    db_gene = crud.get_gene(db, gene_id=gene_id)
    if db_gene is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"유전자 ID {gene_id}를 찾을 수 없습니다."
        )
    
    # 심볼이 변경되는 경우 중복 체크는 crud.update_gene 내에서 처리하거나,
    # 여기에서 미리 체크할 수 있습니다. (현재 crud.py에 update_gene 함수가 정의되어 있지 않으므로, crud.py의 update_gene 로직을 따른다고 가정합니다.)
    # (주의: 이전 CRUD 통합 단계에서 update/delete 함수를 추가했습니다.)
    
    try:
        return crud.update_gene(db, db_gene=db_gene, gene_in=gene_in)
    except ValueError as e:
        # 중복 심볼 오류 처리 (IntegrityError에서 발생)
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


# ======================================================================
# 5. 유전자 삭제 (Delete)
# ======================================================================
@router.delete("/{gene_id}", status_code=status.HTTP_200_OK)
def delete_existing_gene(gene_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 유전자를 삭제합니다. (연관된 변이도 함께 삭제됩니다)"""
    db_gene = crud.get_gene(db, gene_id=gene_id)
    if db_gene is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"유전자 ID {gene_id}를 찾을 수 없습니다."
        )
    
    crud.delete_gene(db=db, db_gene=db_gene)
    return {"message": f"유전자 ID {gene_id} 및 연관 데이터가 성공적으로 삭제되었습니다."}