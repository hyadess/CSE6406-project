"""The explicit inference grid used by the experiment."""

from .gtr import MODEL as GTR
from .gtr_g4 import MODEL as GTR_G4
from .hky import MODEL as HKY
from .hky_g4 import MODEL as HKY_G4
from .jc import MODEL as JC

ANALYSIS_MODELS = (GTR_G4, GTR, HKY_G4, HKY, JC)
MODEL_BY_NAME = {model.name: model for model in ANALYSIS_MODELS}
