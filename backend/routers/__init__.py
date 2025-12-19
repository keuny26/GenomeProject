# Project/routers/__init__.py

# 패키지 외부에서 routers.gene, routers.variant 등으로 접근 가능하도록 모듈을 노출합니다.

from . import document
from . import gene
from . import variant
from . import disease
from . import association

# 또는 더 간단하게, __all__을 사용하여 노출할 수도 있습니다.
# __all__ = ["document", "gene", "variant", "disease", "association"]