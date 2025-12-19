# Project/services/analysis.py

from sqlalchemy.orm import Session
from typing import Dict, Any, List

# 로컬 모듈 임포트
import crud, schemas
from models import DocumentModel

def analyze_document_and_save_graph(db: Session, document_id: int) -> Dict[str, Any]:
    """
    특정 문서의 내용을 AI가 분석하고, 추출된 엔티티와 관계를 DB에 저장합니다.
    (실제 구현에서는 외부 NLP 모델을 호출하는 로직이 들어갑니다.)
    """
    db_document = crud.get_document(db, document_id=document_id)
    if not db_document:
        raise ValueError(f"Document ID {document_id}를 찾을 수 없습니다.")

    # 1. 문서 내용 (db_document.content)을 AI 분석 모델에 전달
    # --- DUMMY AI ANALYSIS RESULT ---
    # 실제 분석에서는 텍스트에서 유전자, 변이, 질병 등을 찾아냅니다.
    
    # 예시 데이터: BRCA1 유전자의 변이 1개와 폐암(Lung Cancer)의 관계 추출
    extracted_data = {
        "genes": [
            {"symbol": "BRCA1", "description": "Breast cancer type 1 susceptibility protein"},
        ],
        "diseases": [
            {"name": "Lung Cancer", "description": "Malignant lung tumor"},
        ],
        "variants": [
            {"gene_symbol": "BRCA1", "name": "V32M", "effect": "Missense"},
        ],
        "associations": [
            {"variant_name": "V32M", "disease_name": "Lung Cancer", "score": 0.85, "summary": "V32M variant in BRCA1 is associated with increased risk of Lung Cancer."}
        ]
    }
    # ------------------------------

    saved_entities = {
        "genes": [],
        "variants": [],
        "diseases": [],
        "associations": []
    }
    
    # 2. 추출된 엔티티를 DB에 저장 (CRUD 함수 활용)
    
    # 2.1 Diseases 저장 (중복 시 기존 엔티티 반환)
    disease_map = {}
    for d_data in extracted_data["diseases"]:
        db_disease = crud.get_or_create_disease(db, name=d_data["name"], defaults={"description": d_data["description"]})
        saved_entities["diseases"].append(schemas.Disease.model_validate(db_disease).model_dump())
        disease_map[d_data["name"]] = db_disease.id
        
    # 2.2 Genes 저장 (중복 시 기존 엔티티 반환)
    gene_map = {}
    for g_data in extracted_data["genes"]:
        db_gene = crud.get_or_create_gene(db, symbol=g_data["symbol"], defaults={"description": g_data["description"]})
        saved_entities["genes"].append(schemas.Gene.model_validate(db_gene).model_dump())
        gene_map[g_data["symbol"]] = db_gene.id

    # 2.3 Variants 저장 (Gene ID 필요)
    variant_map = {}
    for v_data in extracted_data["variants"]:
        gene_id = gene_map.get(v_data["gene_symbol"])
        if gene_id is None:
            # 유전자 엔티티가 누락되었거나 찾을 수 없으면 건너뜀
            continue 
            
        # 변이 이름과 유전자 ID를 기반으로 변이 생성/조회
        db_variant = crud.get_or_create_variant(
            db, 
            gene_id=gene_id,
            name=v_data["name"], 
            defaults={"effect": v_data["effect"], "document_id": document_id}
        )
        saved_entities["variants"].append(schemas.Variant.model_validate(db_variant).model_dump())
        variant_map[v_data["name"]] = db_variant.id

    # 2.4 Associations 저장 (Variant ID, Disease ID 필요)
    for a_data in extracted_data["associations"]:
        variant_id = variant_map.get(a_data["variant_name"])
        disease_id = disease_map.get(a_data["disease_name"])

        if variant_id is None or disease_id is None:
            # 필수 엔티티가 없으면 건너뜀
            continue

        # 연관성 생성 (이미 존재하는 쌍인지 확인)
        assoc_in = schemas.AssociationCreate(
            variant_id=variant_id, 
            disease_id=disease_id, 
            score=a_data["score"],
            summary=a_data["summary"],
            document_id=document_id # 연관성도 문서와 연결
        )
        db_assoc = crud.get_association_by_pair(db, variant_id=variant_id, disease_id=disease_id)
        
        if not db_assoc:
            db_assoc = crud.create_association(db, association=assoc_in)
        
        saved_entities["associations"].append(schemas.Association.model_validate(db_assoc).model_dump())

    # 3. 문서 분석 상태 업데이트 (선택 사항)
    # db_document.is_analyzed = True
    # db.commit()
    
    return saved_entities