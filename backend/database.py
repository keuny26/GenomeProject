# Project/database.py

import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session

# models.py에서 정의한 Base 클래스를 정확하게 임포트합니다.
from models import Base 

# 환경 변수 로드 시 기본값 설정 (None 방지)
# 🚨🚨🚨 DB_HOST와 DB_DATABASE에 기본값 설정 (중요) 🚨🚨🚨
DB_USER = os.getenv("MYSQL_USER", "root")
DB_PASSWORD = os.getenv("MYSQL_PASSWORD", "")
DB_HOST = os.getenv("MYSQL_HOST", "localhost")  # <--- None이 되지 않도록 기본값 설정
DB_NAME = os.getenv("MYSQL_DATABASE", "genomics_db") # <--- None이 되지 않도록 기본값 설정

# MySQL 데이터베이스 URL 구성
SQLALCHEMY_DATABASE_URL = (
    f"mysql+mysqlconnector://{DB_USER}:{DB_PASSWORD}@{DB_HOST}/{DB_NAME}"
)

# SQLAlchemy 엔진 생성
engine = create_engine(
    SQLALCHEMY_DATABASE_URL, 
    pool_recycle=3600, # MySQL 연결이 닫히는 것을 방지
    # debug print 확인을 위해 echo=True를 임시로 설정할 수 있습니다.
    # echo=True 
)

# SessionLocal (데이터베이스 세션 클래스) 생성
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Dependency (세션 생성 및 닫기) 함수
def get_db():
    """요청 시 데이터베이스 세션을 생성하고, 요청 완료 후 세션을 닫습니다."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()