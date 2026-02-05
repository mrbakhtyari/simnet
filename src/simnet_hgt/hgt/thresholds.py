"""Threshold computation and persistence for HGT analysis.

This module provides functions for saving and loading HGT thresholds,
as well as generating analysis reports.
"""

import json
import logging
from pathlib import Path

from .models import (
    KS_EXCELLENT_THRESHOLD,
    KS_FAIR_THRESHOLD,
    KS_GOOD_THRESHOLD,
    BootstrapResult,
    DistributionFit,
    HgtStatistics,
)

logger = logging.getLogger(__name__)


def save_thresholds(bootstrap_result: BootstrapResult, output_path: Path) -> None:
    """Save bootstrap results and determined thresholds to a JSON file.

    Args:
        bootstrap_result: Bootstrap analysis results containing CI bounds.
        output_path: Path to save the JSON file.
    """
    threshold_data = {
        "hgt_95th_percentile_point_estimate": bootstrap_result.point_estimate,
        "hgt_95ci_lower_bound": bootstrap_result.lower_ci,
        "hgt_95ci_upper_bound": bootstrap_result.upper_ci,
        "description": (
            "95% confidence interval for the 95th percentile based on 10,000 bootstrap samples"
        ),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(threshold_data, f, indent=4)

    logger.info("Thresholds saved to %s", output_path)


def load_thresholds(file_path: Path) -> dict:
    """Load previously calculated thresholds from a JSON file.

    Args:
        file_path: Path to the thresholds JSON file.

    Returns:
        Dictionary containing threshold values.

    Raises:
        FileNotFoundError: If the threshold file does not exist.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"No threshold file found at {file_path}")

    with file_path.open(encoding="utf-8") as f:
        return json.load(f)


def _get_fit_quality(ks_statistic: float) -> str:
    """Determine fit quality category based on KS statistic.

    Args:
        ks_statistic: Kolmogorov-Smirnov test statistic.

    Returns:
        Quality category string.
    """
    if ks_statistic < KS_EXCELLENT_THRESHOLD:
        return "Excellent"
    if ks_statistic < KS_GOOD_THRESHOLD:
        return "Good"
    if ks_statistic < KS_FAIR_THRESHOLD:
        return "Fair"
    return "Poor"


def generate_report(  # noqa: PLR0915
    stats: HgtStatistics,
    fitted_dists: list[DistributionFit],
    bootstrap_result: BootstrapResult,
    output_path: Path,
) -> None:
    """Generate a detailed text report with all analysis results.

    Args:
        stats: Descriptive statistics for the HGT scores.
        fitted_dists: List of fitted distribution results.
        bootstrap_result: Bootstrap analysis results.
        output_path: Path to save the report file.
    """
    point_est = bootstrap_result.point_estimate
    lower_ci = bootstrap_result.lower_ci
    upper_ci = bootstrap_result.upper_ci

    with output_path.open("w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("HGT SCORE ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Section 1: Descriptive Statistics
        f.write("1. DESCRIPTIVE STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"  {'count':20s}: {stats.count:15d}\n")
        f.write(f"  {'mean':20s}: {stats.mean:15.6f}\n")
        f.write(f"  {'std':20s}: {stats.std:15.6f}\n")
        f.write(f"  {'median':20s}: {stats.median:15.6f}\n")
        f.write(f"  {'min':20s}: {stats.min:15.6f}\n")
        f.write(f"  {'max':20s}: {stats.max:15.6f}\n")
        f.write(f"  {'q25':20s}: {stats.q25:15.6f}\n")
        f.write(f"  {'q75':20s}: {stats.q75:15.6f}\n")
        f.write(f"  {'q95':20s}: {stats.q95:15.6f}\n")
        f.write(f"  {'q99':20s}: {stats.q99:15.6f}\n")
        f.write(f"  {'iqr':20s}: {stats.iqr:15.6f}\n")
        f.write(f"  {'skewness':20s}: {stats.skewness:15.6f}\n")
        f.write(f"  {'kurtosis':20s}: {stats.kurtosis:15.6f}\n")

        # Section 2: Distribution Fitting
        f.write("\n2. THEORETICAL DISTRIBUTION FITTING\n")
        f.write("-" * 80 + "\n")
        f.write("Kolmogorov-Smirnov Test Results:\n")
        f.write(f"  {'Distribution':<15} {'KS Statistic':<15} {'p-value':<15} {'Fit Quality'}\n")
        f.write("  " + "-" * 70 + "\n")

        for dist in sorted(fitted_dists, key=lambda x: x.ks_statistic):
            fit_quality = _get_fit_quality(dist.ks_statistic)
            f.write(
                f"  {dist.name:<15} {dist.ks_statistic:<15.6f} "
                f"{dist.ks_pvalue:<15.6e} {fit_quality}\n"
            )

        best_fit = min(fitted_dists, key=lambda x: x.ks_statistic)
        f.write(f"\n  Best fitting distribution: {best_fit.name.upper()}\n")
        f.write(f"  Parameters: {best_fit.params}\n")

        # Section 3: Bootstrap Analysis
        f.write("\n3. BOOTSTRAP ANALYSIS (10,000 samples)\n")
        f.write("-" * 80 + "\n")
        f.write("  Target quantile: 95th percentile\n")
        f.write(f"  Point estimate: {point_est:.6f}\n")
        f.write(f"  95% Confidence Interval: [{lower_ci:.6f}, {upper_ci:.6f}]\n")
        f.write(f"  CI Width: {upper_ci - lower_ci:.6f}\n")

        # Section 4: Recommended Threshold
        f.write("\n4. RECOMMENDED HGT THRESHOLD\n")
        f.write("-" * 80 + "\n")
        f.write(f"  Threshold Value: {upper_ci:.6f}\n")
        f.write("  Rationale: Upper bound of 95% CI for 95th percentile\n")
        f.write("  \n")
        f.write("  This threshold ensures that pairs with HGT_score above this value\n")
        f.write("  are significantly higher than expected with 95% confidence.\n")
        f.write("  \n")
        f.write(f"  Pairs to consider for HGT: HGT_score > {upper_ci:.6f}\n")

        # Section 5: Interpretation Guide
        f.write("\n5. INTERPRETATION GUIDE\n")
        f.write("-" * 80 + "\n")
        f.write(f"  • Scores > {upper_ci:.6f}: Likely HGT candidates (high confidence)\n")
        f.write(
            f"  • Scores {point_est:.6f} - {upper_ci:.6f}: Potential HGT (moderate confidence)\n"
        )
        f.write(f"  • Scores < {point_est:.6f}: Less likely HGT events\n")

        f.write("\n" + "=" * 80 + "\n")

    logger.info("Saved detailed report to %s", output_path)
