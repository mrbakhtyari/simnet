"""Data classes, type definitions, and constants for HGT analysis."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np
from numpy.typing import NDArray

# Type aliases for domain clarity
ProteinId = str
GenomeAccession = str
ChgId = str
TaxId = str

# Mapping types
PhylogenyDistances = dict[tuple[GenomeAccession, GenomeAccession], int]
ProteinToGenome = dict[ProteinId, GenomeAccession]
ProteinToChg = dict[ProteinId, ChgId]
GenomeToTaxId = dict[GenomeAccession, TaxId]

# Statistics types
SimilarityStats = dict[int, tuple[float, float]]  # N_int -> (mean, std)

# Constants for distribution fit quality thresholds
KS_EXCELLENT_THRESHOLD: float = 0.05
KS_GOOD_THRESHOLD: float = 0.1
KS_FAIR_THRESHOLD: float = 0.15


class ScoringMethod(str, Enum):
    """HGT scoring method to use.

    CLASSIC: Original formula (N_int * s_AB) / (N_int_mean * s_mean)
    ZSCORE: Z-score based on deviation from expected similarity at phylogenetic distance
    MODIFIED: (s_AB / E[s|N_int]) * log(N_int + 1)
    COMBINED: Z-score weighted by normalized phylogenetic distance
    """

    CLASSIC = "classic"
    ZSCORE = "zscore"
    MODIFIED = "modified"
    COMBINED = "combined"
    EXPONENTIAL = "exponential"
    RANKING = "ranking"
    RANKING_TH = "ranking_th"


@dataclass(frozen=True)
class HgtParams:
    """Parameters for HGT score computation.

    Attributes:
        n_int_mean: Mean phylogenetic distance for normalization.
        hits_mean: Mean protein similarity (pident) for normalization.
        min_phylogenetic_distance: Minimum phylogenetic distance required
            for a hit to be considered. Defaults to 2.
        scoring_method: Which scoring method to use. Defaults to ZSCORE.
    """

    n_int_mean: float
    hits_mean: float
    min_phylogenetic_distance: int = 2
    scoring_method: ScoringMethod = ScoringMethod.ZSCORE
    rank_threshold: float = 0.05


@dataclass(frozen=True)
class BootstrapResult:
    """Result from bootstrap quantile analysis.

    Attributes:
        point_estimate: The observed quantile value from the original data.
        lower_ci: Lower bound of the confidence interval.
        upper_ci: Upper bound of the confidence interval.
        distribution: Array of bootstrap quantile samples.
    """

    point_estimate: float
    lower_ci: float
    upper_ci: float
    distribution: NDArray[np.floating]


@dataclass(frozen=True)
class DistributionFit:
    """Result from fitting a theoretical distribution.

    Attributes:
        name: Name of the distribution (e.g., 'normal', 'lognormal').
        params: Fitted distribution parameters.
        ks_statistic: Kolmogorov-Smirnov test statistic.
        ks_pvalue: Kolmogorov-Smirnov test p-value.
        distribution: The scipy distribution object.
    """

    name: str
    params: tuple
    ks_statistic: float
    ks_pvalue: float
    distribution: Any  # scipy.stats.rv_continuous or rv_frozen


@dataclass
class ChgProperties:
    """Properties of a Cluster of Homologous Genes (CHG).

    Attributes:
        chg_id: Unique identifier for the CHG.
        genome_set: Set of genome accessions containing proteins in this CHG.
        max_nint: Maximum phylogenetic distance between any pair of genomes in the CHG.
    """

    chg_id: ChgId
    genome_set: set[GenomeAccession] = field(default_factory=set)
    max_nint: int = -1


@dataclass
class HgtStatistics:
    """Descriptive statistics for HGT scores.

    Attributes:
        count: Number of scores.
        mean: Mean score.
        std: Standard deviation.
        median: Median score.
        min: Minimum score.
        max: Maximum score.
        q25: 25th percentile.
        q75: 75th percentile.
        q95: 95th percentile.
        q99: 99th percentile.
        iqr: Interquartile range.
        skewness: Skewness of the distribution.
        kurtosis: Kurtosis of the distribution.
    """

    count: int
    mean: float
    std: float
    median: float
    min: float
    max: float
    q25: float
    q75: float
    q95: float
    q99: float
    iqr: float
    skewness: float
    kurtosis: float
