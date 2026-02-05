"""HGT (Horizontal Gene Transfer) analysis module.

This module provides functions for computing and analyzing HGT scores
from protein similarity data and phylogenetic distances.
"""

from .analyze_scores import analyze_hgt_scores
from .models import HgtParams, ScoringMethod
from .postprocessing import restore_original_names
from .score import compute_hgt_scores, compute_similarity_mean

__all__ = [
    "HgtParams",
    "ScoringMethod",
    "analyze_hgt_scores",
    "compute_hgt_scores",
    "compute_similarity_mean",
    "restore_original_names",
]
