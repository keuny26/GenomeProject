# Project/routers/association.py

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import List, Optional

# 로컬 모듈 임포트
import crud, schemas
from database import get_db
from models import AssociationModel

# 라우터 인스턴스 생성
router = APIRouter(
    prefix="/associations",
    tags=["Associations (연관성)"],
)

# ======================================================================
# 1. 연관성 목록 조회 (Read All)
# ======================================================================
@router.get("/", response_model=List[schemas.Association])
def read_associations(
    variant_id: Optional[int] = Query(None, description="특정 변이 ID로 필터링"),
    disease_id: Optional[int] = Query(None, description="특정 질병 ID로 필터링"),
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """
    모든 연관성 목록을 조회하거나, 변이/질병 ID로 필터링하여 조회합니다.
    """
    associations = crud.get_associations(
        db, 
        variant_id=variant_id, 
        disease_id=disease_id, 
        skip=skip, 
        limit=limit
    )
    return associations

# ======================================================================
# 2. 연관성 생성 (Create)
# ======================================================================
@router.post("/", response_model=schemas.Association, status_code=status.HTTP_201_CREATED)
def create_new_association(association: schemas.AssociationCreate, db: Session = Depends(get_db)):
    """
    새로운 연관성을 생성합니다. (Variant ID와 Disease ID 쌍은 중복될 수 없습니다.)
    """
    # 중복 체크 (DB IntegrityError를 이용해도 되지만, 명시적 체크로 오류 메시지 명확화)
    existing_association = crud.get_association_by_pair(
        db, 
        variant_id=association.variant_id, 
        disease_id=association.disease_id
    )
    if existing_association:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Association between Variant ID {association.variant_id} and Disease ID {association.disease_id} already exists."
        )
        
    try:
        # crud.create_association 함수 내에서 FK 유효성 검사 (Variant/Disease ID 존재 여부)를 수행합니다.
        return crud.create_association(db=db, association=association)
    except ValueError as e:
        # FK 오류 또는 기타 DB 오류 처리
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, 
            detail=str(e)
        )


# ======================================================================
# 3. 특정 연관성 조회 (Read One)
# ======================================================================
@router.get("/{association_id}", response_model=schemas.Association)
def read_association(association_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 연관성을 조회합니다."""
    db_association = crud.get_association(db, association_id=association_id)
    if db_association is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"연관성 ID {association_id}를 찾을 수 없습니다."
        )
    return db_association

# ======================================================================
# 4. 연관성 수정 (Update)
# ======================================================================
@router.put("/{association_id}", response_model=schemas.Association)
def update_existing_association(association_id: int, association_in: schemas.AssociationCreate, db: Session = Depends(get_db)):
    """특정 ID를 가진 연관성 정보를 수정합니다."""
    db_association = crud.get_association(db, association_id=association_id)
    if db_association is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"연관성 ID {association_id}를 찾을 수 없습니다."
        )
    
    try:
        # crud.update_association는 이미 crud.py에 정의되어 있습니다.
        return crud.update_association(db, db_association=db_association, association_in=association_in)
    except ValueError as e:
        # FK 오류 또는 중복 조합 오류 처리
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )

# ======================================================================
# 5. 연관성 삭제 (Delete)
# ======================================================================
@router.delete("/{association_id}", status_code=status.HTTP_200_OK)
def delete_existing_association(association_id: int, db: Session = Depends(get_db)):
    """특정 ID를 가진 연관성을 삭제합니다."""
    db_association = crud.get_association(db, association_id=association_id)
    if db_association is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, 
            detail=f"연관성 ID {association_id}를 찾을 수 없습니다."
        )
    
    crud.delete_association(db=db, db_association=db_association)
    return {"message": f"연관성 ID {association_id}가 성공적으로 삭제되었습니다."}