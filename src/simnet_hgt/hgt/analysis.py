"""Statistical analysis functions for HGT scores.

This module provides functions for computing descriptive statistics,
fitting theoretical distributions, and bootstrap analysis for HGT scores.
"""

import logging

import numpy as np
from numpy.typing import NDArray
from scipy import stats
from scipy.stats import expon, gamma, lognorm, norm
from tqdm import tqdm

from .models import BootstrapResult, DistributionFit, HgtStatistics

logger = logging.getLogger(__name__)


def compute_basic_statistics(scores: NDArray[np.floating]) -> HgtStatistics:
    """Calculate descriptive statistics for HGT scores.

    Args:
        scores: Array of HGT score values.

    Returns:
        HgtStatistics dataclass containing all computed statistics.
    """
    return HgtStatistics(
        count=len(scores),
        mean=float(np.mean(scores)),
        std=float(np.std(scores)),
        median=float(np.median(scores)),
        min=float(np.min(scores)),
        max=float(np.max(scores)),
        q25=float(np.percentile(scores, 25)),
        q75=float(np.percentile(scores, 75)),
        q95=float(np.percentile(scores, 95)),
        q99=float(np.percentile(scores, 99)),
        iqr=float(np.percentile(scores, 75) - np.percentile(scores, 25)),
        skewness=float(stats.skew(scores)),
        kurtosis=float(stats.kurtosis(scores)),
    )


def fit_distributions(scores: NDArray[np.floating]) -> list[DistributionFit]:
    """Fit theoretical distributions to the HGT score data.

    Fits normal, log-normal, gamma, and exponential distributions
    and evaluates fit quality using the Kolmogorov-Smirnov test.

    Args:
        scores: Array of HGT score values.

    Returns:
        List of DistributionFit dataclasses sorted by KS statistic (best first).
    """
    logger.info("Fitting theoretical distributions...")
    fitted_dists: list[DistributionFit] = []

    # Normal distribution
    mu, sigma = norm.fit(scores)
    ks_stat, ks_pval = stats.kstest(scores, "norm", args=(mu, sigma))
    fitted_dists.append(
        DistributionFit(
            name="normal",
            params=(mu, sigma),
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            distribution=norm,
        )
    )

    # Log-normal distribution (requires positive values)
    if np.all(scores > 0):
        shape, loc, scale = lognorm.fit(scores, floc=0)
        ks_stat, ks_pval = stats.kstest(scores, "lognorm", args=(shape, loc, scale))
        fitted_dists.append(
            DistributionFit(
                name="lognormal",
                params=(shape, loc, scale),
                ks_statistic=ks_stat,
                ks_pvalue=ks_pval,
                distribution=lognorm,
            )
        )

    # Gamma distribution (requires positive values)
    if np.all(scores > 0):
        shape, loc, scale = gamma.fit(scores)
        ks_stat, ks_pval = stats.kstest(scores, "gamma", args=(shape, loc, scale))
        fitted_dists.append(
            DistributionFit(
                name="gamma",
                params=(shape, loc, scale),
                ks_statistic=ks_stat,
                ks_pvalue=ks_pval,
                distribution=gamma,
            )
        )

    # Exponential distribution (requires positive values)
    if np.all(scores > 0):
        loc, scale = expon.fit(scores)
        ks_stat, ks_pval = stats.kstest(scores, "expon", args=(loc, scale))
        fitted_dists.append(
            DistributionFit(
                name="exponential",
                params=(loc, scale),
                ks_statistic=ks_stat,
                ks_pvalue=ks_pval,
                distribution=expon,
            )
        )

    # Sort by KS statistic (lower is better)
    fitted_dists.sort(key=lambda x: x.ks_statistic)

    for dist in fitted_dists:
        logger.info(
            "  %s: KS statistic = %.6f, p-value = %.6e",
            dist.name,
            dist.ks_statistic,
            dist.ks_pvalue,
        )

    return fitted_dists


def bootstrap_quantile(
    scores: NDArray[np.floating],
    quantile: float = 0.95,
    n_bootstrap: int = 10000,
    confidence_level: float = 0.95,
) -> BootstrapResult:
    """Perform bootstrap analysis to determine confidence interval for a quantile.

    Uses resampling with replacement to estimate the sampling distribution
    of a quantile and compute confidence intervals.

    Args:
        scores: Array of HGT score values.
        quantile: Target quantile to estimate (default: 0.95 for 95th percentile).
        n_bootstrap: Number of bootstrap samples (default: 10,000).
        confidence_level: Confidence level for the interval (default: 0.95).

    Returns:
        BootstrapResult containing point estimate, confidence bounds, and
        the full bootstrap distribution.
    """
    logger.info("Performing bootstrap analysis with %d samples...", n_bootstrap)

    bootstrap_quantiles: list[float] = []
    n = len(scores)
    show_progress = logger.isEnabledFor(logging.INFO)

    for bootstrap_sample in tqdm(
        (np.random.choice(scores, size=n, replace=True) for _ in range(n_bootstrap)),
        total=n_bootstrap,
        desc="Bootstrap sampling",
        disable=not show_progress,
    ):
        bootstrap_quantiles.append(float(np.percentile(bootstrap_sample, quantile * 100)))

    bootstrap_array = np.array(bootstrap_quantiles)

    # Calculate confidence interval
    alpha = 1 - confidence_level
    lower_bound = float(np.percentile(bootstrap_array, alpha / 2 * 100))
    upper_bound = float(np.percentile(bootstrap_array, (1 - alpha / 2) * 100))
    point_estimate = float(np.percentile(scores, quantile * 100))

    logger.info(
        "Bootstrap results for %.0f%% percentile: point=%.6f, %.0f%% CI=[%.6f, %.6f]",
        quantile * 100,
        point_estimate,
        confidence_level * 100,
        lower_bound,
        upper_bound,
    )

    return BootstrapResult(
        point_estimate=point_estimate,
        lower_ci=lower_bound,
        upper_ci=upper_bound,
        distribution=bootstrap_array,
    )


def get_best_fit(distributions: list[DistributionFit]) -> DistributionFit:
    """Get the best fitting distribution from a list.

    Args:
        distributions: List of fitted distributions.

    Returns:
        The distribution with the lowest KS statistic.

    Raises:
        ValueError: If the list is empty.
    """
    if not distributions:
        raise ValueError("No distributions provided")
    return min(distributions, key=lambda x: x.ks_statistic)
