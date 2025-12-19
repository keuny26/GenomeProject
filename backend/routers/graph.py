from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import select
from typing import Dict, List, Any

# 로컬 모듈 임포트
import crud
from database import get_db
from models import GeneModel, VariantModel, AssociationModel, DiseaseModel
import schemas 

router = APIRouter(
    prefix="/graph",
    tags=["Graph Lookups (그래프 조회)"],
)

# --- (기존 1번 유전자 조회 / 2번 변이 조회 코드는 그대로 유지) ---

# [1. 특정 유전자 ID 기반 조회 코드...]
# [2. 특정 변이 ID 기반 조회 코드...]

# ======================================================================
# 3. [추가] 지식 그래프 시각화용 전체 데이터 조회
# ======================================================================

@router.get("/data", response_model=Dict[str, Any])
def get_full_graph_data(db: Session = Depends(get_db)):
    """
    프론트엔드 Force-Graph 시각화를 위해 모든 Gene, Disease 노드와 
    그 사이의 Association 링크를 반환합니다.
    """
    # 1. 모든 노드(유전자, 질병) 및 관계 로드
    genes = db.query(GeneModel).all()
    diseases = db.query(DiseaseModel).all()
    # 관계 조회를 위해 Variant 정보까지 로드
    associations = db.query(AssociationModel).options(
        selectinload(AssociationModel.variant)
    ).all()

    nodes = []
    links = []

    # 2. 유전자 노드 생성
    for g in genes:
        nodes.append({
            "id": f"gene_{g.id}",
            "label": f"🧬 {g.symbol}",
            "type": "gene",
            "color": "#4285F4"  # 파란색
        })

    # 3. 질병 노드 생성
    for d in diseases:
        nodes.append({
            "id": f"disease_{d.id}",
            "label": f"🏥 {d.name}",
            "type": "disease",
            "color": "#EA4335"  # 빨간색
        })

    # 4. 연결선(Link) 생성 (Gene <-> Disease)
    # Association은 Variant를 통해 Gene과 연결됩니다.
    for assoc in associations:
        if assoc.variant and assoc.variant.gene_id:
            links.append({
                "source": f"gene_{assoc.variant.gene_id}",
                "target": f"disease_{assoc.disease_id}",
                "value": assoc.p_value or 1,
                "label": assoc.variant.name  # 선 위에 변이 이름 표시
            })

    return {"nodes": nodes, "links": links}