from sqlalchemy.orm import Session
from sqlalchemy import select
from typing import List, Optional
import models 
import schemas 
from sqlalchemy.exc import IntegrityError

# ======================================================================
# 0. Document (문서) CRUD
# ======================================================================
def get_document(db: Session, document_id: int) -> Optional[models.DocumentModel]:
    return db.scalar(select(models.DocumentModel).where(models.DocumentModel.id == document_id))
    
def get_documents(db: Session, skip: int = 0, limit: int = 100) -> List[models.DocumentModel]:
    return db.scalars(select(models.DocumentModel).offset(skip).limit(limit)).all()

def create_document(db: Session, file_name: str, content: str) -> models.DocumentModel:
    db_document = models.DocumentModel(name=file_name, content=content) 
    db.add(db_document)
    db.flush() 
    return db_document

# ======================================================================
# 1. Gene (유전자) CRUD
# ======================================================================
def get_gene(db: Session, gene_id: int) -> Optional[models.GeneModel]:
    return db.scalar(select(models.GeneModel).where(models.GeneModel.id == gene_id))

def get_genes(db: Session, skip: int = 0, limit: int = 100) -> List[models.GeneModel]:
    return db.scalars(select(models.GeneModel).offset(skip).limit(limit)).all()

def get_gene_by_symbol(db: Session, symbol: str) -> Optional[models.GeneModel]:
    return db.scalar(select(models.GeneModel).where(models.GeneModel.symbol == symbol).limit(1))

def get_or_create_gene(db: Session, symbol: str, name: Optional[str] = None, location: Optional[str] = None) -> models.GeneModel:
    db_gene = get_gene_by_symbol(db, symbol)
    if db_gene: return db_gene
    
    db_gene = models.GeneModel(symbol=symbol, name=name, location=location)
    db.add(db_gene)
    db.flush() 
    return db_gene

def create_gene(db: Session, gene: schemas.GeneCreate) -> models.GeneModel:
    db_gene = models.GeneModel(symbol=gene.symbol, name=gene.name, location=gene.location)
    db.add(db_gene)
    try:
        db.flush()
    except IntegrityError:
        db.rollback()
        raise ValueError(f"Gene '{gene.symbol}' already exists.")
    return db_gene

# ======================================================================
# 2. Variant (변이) CRUD
# ======================================================================
def get_variant(db: Session, variant_id: int) -> Optional[models.VariantModel]:
    return db.scalar(select(models.VariantModel).where(models.VariantModel.id == variant_id))

def get_variants(db: Session, gene_id: Optional[int] = None, skip: int = 0, limit: int = 100) -> List[models.VariantModel]:
    query = select(models.VariantModel)
    if gene_id is not None:
        query = query.where(models.VariantModel.gene_id == gene_id)
    return db.scalars(query.offset(skip).limit(limit)).all()

def create_variant(db: Session, variant: schemas.VariantCreate) -> models.VariantModel:
    db_variant = models.VariantModel(
        gene_id=variant.gene_id, name=variant.name, 
        variant_type=variant.variant_type, document_id=variant.document_id
    )
    db.add(db_variant)
    try:
        db.flush()
    except IntegrityError:
        db.rollback()
        raise ValueError(f"Variant '{variant.name}' already exists.")
    return db_variant

# ======================================================================
# 3. Disease (질병) CRUD
# ======================================================================
def get_disease_by_name(db: Session, name: str) -> Optional[models.DiseaseModel]:
    return db.scalar(select(models.DiseaseModel).where(models.DiseaseModel.name == name).limit(1))

def get_or_create_disease(db: Session, name: str, description: Optional[str] = None) -> models.DiseaseModel:
    db_disease = get_disease_by_name(db, name)
    if db_disease: return db_disease

    db_disease = models.DiseaseModel(name=name, description=description)
    db.add(db_disease)
    try:
        db.flush()
        return db_disease
    except IntegrityError:
        db.rollback()
        return get_disease_by_name(db, name)

# ======================================================================
# 4. Association (연관성) CRUD
# ======================================================================
def create_association(db: Session, association: schemas.AssociationCreate) -> models.AssociationModel:
    db_association = models.AssociationModel(
        variant_id=association.variant_id, disease_id=association.disease_id,
        p_value=association.p_value, odds_ratio=association.odds_ratio,
        reference=association.reference, document_id=association.document_id
    )
    db.add(db_association)
    try:
        db.flush()
    except IntegrityError:
        db.rollback()
        raise ValueError("Association already exists.")
    return db_association

def delete_gene(db: Session, db_gene: models.GeneModel):
    db.delete(db_gene)
    db.commit()

def delete_variant(db: Session, db_variant: models.VariantModel):
    db.delete(db_variant)
    db.commit()

def delete_disease(db: Session, db_disease: models.DiseaseModel):
    db.delete(db_disease)
    db.commit()