"""HGT scoring strategies.

This module provides scoring functions for computing HGT scores using
different methodologies. Each function follows a common signature pattern.
"""

import math
from collections.abc import Callable

from .models import ScoringMethod

# Type alias for similarity statistics per phylogenetic distance
SimilarityStats = dict[int, tuple[float, float]]  # N_int -> (mean, std)

# Type alias for scorer function signature
Scorer = Callable[..., float]


def calculate_classic_score(
    n_int: int,
    similarity: float,
    n_int_mean: float,
    s_mean: float,
    **_kwargs,
) -> float:
    """Calculate classic HGT score: (N_int * similarity) / (N_int_mean * s_mean).

    Args:
        n_int: Phylogenetic distance (internal nodes) between genomes.
        similarity: Protein similarity (percent identity).
        n_int_mean: Mean phylogenetic distance for normalization.
        s_mean: Mean protein similarity for normalization.

    Returns:
        HGT score value. Returns 0.0 if normalization values are invalid.
    """
    if n_int_mean <= 0 or s_mean <= 0:
        return 0.0
    return (n_int * similarity) / (n_int_mean * s_mean)


def calculate_zscore(
    similarity: float,
    n_int: int,
    similarity_stats: SimilarityStats,
    **_kwargs,
) -> float:
    """Calculate z-score based HGT score.

    Measures how many standard deviations the observed similarity
    is above the expected similarity for the given phylogenetic distance.

    Args:
        similarity: Observed protein similarity (percent identity).
        n_int: Phylogenetic distance (internal nodes).
        similarity_stats: Dictionary mapping N_int to (mean, std).

    Returns:
        Z-score value. Higher values indicate potential HGT.
    """
    if n_int not in similarity_stats:
        return 0.0

    expected_mean, expected_std = similarity_stats[n_int]
    if expected_std <= 0:
        return 0.0

    return (similarity - expected_mean) / expected_std


def calculate_modified_score(
    similarity: float,
    n_int: int,
    similarity_stats: SimilarityStats,
    **_kwargs,
) -> float:
    """Calculate modified HGT score: (s_AB / E[s|N_int]) * log(N_int + 1).

    Uses ratio of observed to expected similarity, weighted by log-scaled distance.

    Args:
        similarity: Observed protein similarity (percent identity).
        n_int: Phylogenetic distance (internal nodes).
        similarity_stats: Dictionary mapping N_int to (mean, std).

    Returns:
        Modified HGT score. Higher values indicate potential HGT.
    """
    if n_int not in similarity_stats:
        return 0.0

    expected_mean, _ = similarity_stats[n_int]
    if expected_mean <= 0:
        return 0.0

    return (similarity / expected_mean) * math.log(n_int + 1)


def calculate_combined_score(
    similarity: float,
    n_int: int,
    similarity_stats: SimilarityStats,
    n_int_max: int = 10,
    **_kwargs,
) -> float:
    """Calculate combined HGT score: z_score * (1 + log(N_int) / log(N_int_max)).

    Combines z-score with normalized phylogenetic distance weighting.

    Args:
        similarity: Observed protein similarity (percent identity).
        n_int: Phylogenetic distance (internal nodes).
        similarity_stats: Dictionary mapping N_int to (mean, std).
        n_int_max: Maximum phylogenetic distance for normalization.

    Returns:
        Combined HGT score. Higher values indicate potential HGT.
    """
    z = calculate_zscore(similarity, n_int, similarity_stats)
    if n_int <= 0 or n_int_max <= 1:
        return z

    distance_weight = 1 + math.log(n_int) / math.log(n_int_max)
    return z * distance_weight


def calculate_exponential_score(  # noqa: PLR0913
    similarity: float,
    n_int: int,
    similarity_stats: SimilarityStats,
    taxon_stats: dict[str, tuple[float, float]] | None = None,
    tax_a: str | None = None,
    tax_b: str | None = None,
    n_cdf: dict[int, float] | None = None,
    exponent: float = 1.0,
    **_kwargs,
) -> float:
    """Calculate HGT score using Taxon-Specific Null Model.

    Instead of global N_int buckets, we compare s_AB to the expected similarity
    for these specific organisms (background vertical inheritance).

    Score = Local_Z_Score * Sigmoid_Distance_Weight
    """
    # Calculate Local Z-Score
    z_score = 0.0
    if taxon_stats and tax_a and tax_b and tax_a in taxon_stats and tax_b in taxon_stats:
        mean_a, std_a = taxon_stats[tax_a]
        mean_b, std_b = taxon_stats[tax_b]

        # Expected similarity is average of background for both taxa
        expected_s = (mean_a + mean_b) / 2.0
        # Pooled std (conservative)
        pooled_std = max(std_a, std_b, 5.0)

        z_score = (similarity - expected_s) / pooled_std
    else:
        # Fallback to Global N_int Z-score if stats missing
        z_score = calculate_zscore(similarity, n_int, similarity_stats)

    if z_score < 0:
        return 0.0

    # Distance Weighting: Use Percentile Rank (CDF) of N_int.
    weight = 1.0
    if n_cdf:
        weight = 1.0 + n_cdf.get(n_int, 0.0)

    # Use exponent to control "Deep Preference".
    return z_score * (weight**exponent)


def calculate_ranking_score(
    n_int: int,
    similarity: float,
    rank: float,
    **_kwargs,
) -> float:
    """Calculate Ranking-based HGT score.

    Formula: (N_int^2 * similarity) / sqrt(rank)

    Args:
        n_int: Phylogenetic distance.
        similarity: Protein similarity.
        rank: Rank of the similarity in the local network (one-way).

    Returns:
        HGT score.
    """
    if rank <= 0:
        return 0.0
    return (math.pow(n_int, 2) * similarity) / math.sqrt(rank)


def calculate_ranking_th_score(
    n_int: int,
    similarity: float,
    rank: float,
    **_kwargs,
) -> float:
    """Calculate Threshold-based Ranking HGT score.

    Formula: (N_int * similarity) / rank
    """
    if rank <= 0:
        return 0.0
    if n_int <= 1:
        return 0.0

    return ((n_int + 1) * similarity) / rank


# Scorer registry mapping methods to functions
_SCORER_REGISTRY: dict[ScoringMethod, Scorer] = {
    ScoringMethod.CLASSIC: calculate_classic_score,
    ScoringMethod.ZSCORE: calculate_zscore,
    ScoringMethod.MODIFIED: calculate_modified_score,
    ScoringMethod.COMBINED: calculate_combined_score,
    ScoringMethod.EXPONENTIAL: calculate_exponential_score,
    ScoringMethod.RANKING: calculate_ranking_score,
    ScoringMethod.RANKING_TH: calculate_ranking_th_score,
}


def get_scorer(method: ScoringMethod) -> Scorer:
    """Get the scoring function for a given method.

    Args:
        method: The scoring method to use.

    Returns:
        The scoring function for the given method.

    Raises:
        ValueError: If the method is not recognized.
    """
    if method not in _SCORER_REGISTRY:
        raise ValueError(f"Unknown scoring method: {method}")
    return _SCORER_REGISTRY[method]
