# Project/routers/disease.py

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import List, Optional

# 로컬 모듈 임포트
import crud, schemas
from database import get_db
from models import DiseaseModel

# 라우터 인스턴스 생성
router = APIRouter(
    prefix="/diseases",
    tags=["Diseases (질병)"],
)

# ======================================================================
# 1. 질병 목록 조회 (Read All)
# ======================================================================
@router.get("/", response_model=List[schemas.Disease])
def read_diseases(
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """
    모든 질병 목록을 조회합니다.
    """
    diseases = crud.get_diseases(db, skip=skip, limit=limit)
    return diseases

# ======================================================================
# 2. 질병 생성 (Create)
# ======================================================================
@router.post("/", response_model=schemas.Disease, status_code=status.HTTP_201_CREATED)
def create_new_disease(disease: schemas.DiseaseCreate, db: Session = Depends(get_db)):
    """
    새로운 질병을 생성합니다. (질병 이름은 중복될 수 없습니다.)
    """
    # 중복 체크
    existing_disease = crud.get_disease_by_name(db, name=disease.name)
    if existing_disease:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Disease with name '{disease.name}' already exists."
        )

    # crud.py에 정의된 get_or_create_disease 함수를 활용하여 생성 (이미 중복 체크를 했으므로 생성될 것임)
    try:
        new_disease = crud.get_or_create_disease(db, name=disease.name)
        return new_disease
    except Exception as e:
         raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, 
            detail=f"Disease creation failed due to database error: {str(e)}"
        )


# ======================================================================
# 3. 특정 질병 조회 (Read One)
# ======================================================================
@router.get("/{disease_id}", response_model=schemas.Disease)
def read_disease(disease_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 질병을 조회합니다."""
    db_disease = crud.get_disease(db, disease_id=disease_id)
    if db_disease is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"질병 ID {disease_id}를 찾을 수 없습니다."
        )
    return db_disease

# ======================================================================
# 4. 질병 수정 (Update)
# ======================================================================
@router.put("/{disease_id}", response_model=schemas.Disease)
def update_existing_disease(disease_id: int, disease_in: schemas.DiseaseCreate, db: Session = Depends(get_db)):
    """특정 ID를 가진 질병 정보를 수정합니다."""
    db_disease = crud.get_disease(db, disease_id=disease_id)
    if db_disease is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"질병 ID {disease_id}를 찾을 수 없습니다."
        )
    
    try:
        # crud.update_disease는 이미 crud.py에 정의되어 있습니다.
        return crud.update_disease(db, db_disease=db_disease, disease_in=disease_in)
    except ValueError as e:
        # 중복 이름 오류 처리 (IntegrityError에서 발생)
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )

# ======================================================================
# 5. 질병 삭제 (Delete)
# ======================================================================
@router.delete("/{disease_id}", status_code=status.HTTP_200_OK)
def delete_existing_disease(disease_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 질병을 삭제합니다. (연관된 Association도 함께 삭제됩니다)"""
    db_disease = crud.get_disease(db, disease_id=disease_id)
    if db_disease is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"질병 ID {disease_id}를 찾을 수 없습니다."
        )
    
    crud.delete_disease(db=db, db_disease=db_disease)
    return {"message": f"질병 ID {disease_id} 및 연관 데이터가 성공적으로 삭제되었습니다."}