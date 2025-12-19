# Project/routers/variant.py

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import List, Optional

# 로컬 모듈 임포트
import crud, schemas
from database import get_db
from models import VariantModel

# 라우터 인스턴스 생성
router = APIRouter(
    prefix="/variants",
    tags=["Variants (변이)"],
)

# ======================================================================
# 1. 변이 목록 조회 (Read All)
# ======================================================================
@router.get("/", response_model=List[schemas.Variant])
def read_variants(
    gene_id: Optional[int] = Query(None, description="특정 유전자에 속한 변이만 필터링"),
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """
    모든 변이 목록을 조회하거나, 특정 유전자(gene_id)에 속한 변이만 필터링하여 조회합니다.
    """
    variants = crud.get_variants(db, gene_id=gene_id, skip=skip, limit=limit)
    return variants

# ======================================================================
# 2. 변이 생성 (Create)
# ======================================================================
@router.post("/", response_model=schemas.Variant, status_code=status.HTTP_201_CREATED)
def create_new_variant(variant: schemas.VariantCreate, db: Session = Depends(get_db)):
    """
    새로운 변이를 생성합니다. (변이는 반드시 존재하는 gene_id에 연결되어야 합니다.)
    """
    try:
        return crud.create_variant(db=db, variant=variant)
    except ValueError as e:
        # FK 오류 또는 중복 조합 오류 (UniqueConstraint) 처리
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, 
            detail=str(e)
        )

# ======================================================================
# 3. 특정 변이 조회 (Read One)
# ======================================================================
@router.get("/{variant_id}", response_model=schemas.Variant)
def read_variant(variant_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 변이를 조회합니다."""
    db_variant = crud.get_variant(db, variant_id=variant_id)
    if db_variant is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"변이 ID {variant_id}를 찾을 수 없습니다."
        )
    return db_variant

# ======================================================================
# 4. 변이 수정 (Update)
# ======================================================================
@router.put("/{variant_id}", response_model=schemas.Variant)
def update_existing_variant(variant_id: int, variant_in: schemas.VariantCreate, db: Session = Depends(get_db)):
    """특정 ID를 가진 변이 정보를 수정합니다."""
    db_variant = crud.get_variant(db, variant_id=variant_id)
    if db_variant is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"변이 ID {variant_id}를 찾을 수 없습니다."
        )
    
    try:
        return crud.update_variant(db, db_variant=db_variant, variant_in=variant_in)
    except ValueError as e:
        # FK 오류 또는 중복 조합 오류 처리
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )

# ======================================================================
# 5. 변이 삭제 (Delete)
# ======================================================================
@router.delete("/{variant_id}", status_code=status.HTTP_200_OK)
def delete_existing_variant(variant_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 변이를 삭제합니다. (연관된 Association도 함께 삭제됩니다)"""
    db_variant = crud.get_variant(db, variant_id=variant_id)
    if db_variant is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"변이 ID {variant_id}를 찾을 수 없습니다."
        )
    
    crud.delete_variant(db=db, db_variant=db_variant)
    return {"message": f"변이 ID {variant_id} 및 연관 데이터가 성공적으로 삭제되었습니다."}