# Project/models.py (SQLAlchemy 2.0 Mapped 최종 버전)

from sqlalchemy import ForeignKey, String, Text, Float, UniqueConstraint
from sqlalchemy.orm import relationship, DeclarativeBase, Mapped, mapped_column
from typing import List, Optional

# Base 클래스 정의
class Base(DeclarativeBase):
    pass

# ======================================================================
# 0. Document (문서) 모델
# ======================================================================
class DocumentModel(Base):
    __tablename__ = "documents"

    # 모든 컬럼을 Mapped 타입으로 정의 (SQLAlchemy 2.0 스타일)
    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    
    # ⚠️ name 컬럼의 중복을 허용 (unique=False)
    name: Mapped[str] = mapped_column(String(255), index=True, comment="파일 이름 또는 제목")
    content: Mapped[str] = mapped_column(Text, comment="문서의 전체 텍스트 내용")

    # 1:N 관계 정의
    associations: Mapped[List["AssociationModel"]] = relationship(
        "AssociationModel", 
        back_populates="document", 
        cascade="all, delete-orphan"
    )

# ======================================================================
# 1. Gene (유전자) 모델
# ======================================================================
class GeneModel(Base):
    __tablename__ = "genes"

    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    symbol: Mapped[str] = mapped_column(String(50), unique=True, index=True)
    name: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    location: Mapped[Optional[str]] = mapped_column(String(100), nullable=True)
    
    # 1:N 관계 정의
    variants: Mapped[List["VariantModel"]] = relationship(
        "VariantModel", 
        back_populates="gene", 
        cascade="all, delete-orphan"
    )

# ======================================================================
# 2. Variant (변이) 모델
# ======================================================================
class VariantModel(Base):
    __tablename__ = "variants"

    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    
    # 외래 키 정의
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id"))
    document_id: Mapped[Optional[int]] = mapped_column(ForeignKey("documents.id"), nullable=True) 

    name: Mapped[str] = mapped_column(String(100), index=True)
    variant_type: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)

    # N:1 관계 정의
    gene: Mapped["GeneModel"] = relationship("GeneModel", back_populates="variants")
    
    # 1:N 관계 정의
    associations: Mapped[List["AssociationModel"]] = relationship(
        "AssociationModel", 
        back_populates="variant", 
        cascade="all, delete-orphan"
    )

# ======================================================================
# 3. Disease (질병) 모델
# ======================================================================
class DiseaseModel(Base):
    __tablename__ = "diseases"

    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    name: Mapped[str] = mapped_column(String(255), unique=True, index=True)
    description: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    # 1:N 관계 정의
    associations: Mapped[List["AssociationModel"]] = relationship(
        "AssociationModel", 
        back_populates="disease", 
        cascade="all, delete-orphan"
    )

# ======================================================================
# 4. Association (연관성) 모델
# ======================================================================
class AssociationModel(Base):
    __tablename__ = "associations"

    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    
    # 외래 키 정의
    variant_id: Mapped[int] = mapped_column(ForeignKey("variants.id"))
    disease_id: Mapped[int] = mapped_column(ForeignKey("diseases.id"))
    document_id: Mapped[Optional[int]] = mapped_column(ForeignKey("documents.id"), nullable=True)

    # 연관성 상세 정보 (Optional로 변경하여 Null 허용)
    score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    p_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    odds_ratio: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    reference: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    
    __table_args__ = (
        UniqueConstraint('variant_id', 'disease_id', name='_variant_disease_uc'),
    )

    # N:1 관계 정의
    variant: Mapped["VariantModel"] = relationship("VariantModel", back_populates="associations")
    disease: Mapped["DiseaseModel"] = relationship("DiseaseModel", back_populates="associations")
    document: Mapped["DocumentModel"] = relationship("DocumentModel", back_populates="associations")