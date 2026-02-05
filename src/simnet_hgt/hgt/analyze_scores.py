"""HGT score analysis pipeline.

This module provides the main analysis pipeline that orchestrates
statistical analysis, visualization, and report generation for HGT scores.
"""

import logging
import warnings
from pathlib import Path

import numpy as np

from .analysis import bootstrap_quantile, compute_basic_statistics, fit_distributions
from .loaders import load_hgt_scores
from .thresholds import generate_report, save_thresholds
from .visualization import create_comprehensive_plots

logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")


def analyze_hgt_scores(
    hgt_scores_tsv: Path,
    output_dir: Path,
    save_threshold_json: bool = True,
) -> None:
    """Main analysis pipeline for HGT scores.

    Performs comprehensive analysis including:
    - Descriptive statistics computation
    - Theoretical distribution fitting
    - Bootstrap analysis for threshold determination
    - Visualization generation
    - Report generation
    - High-confidence HGT candidate identification

    Args:
        hgt_scores_tsv: Path to TSV file containing HGT scores.
        output_dir: Directory to save all output files.
        save_threshold_json: Whether to save thresholds to JSON. Defaults to True.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    df = load_hgt_scores(hgt_scores_tsv)

    # Extract HGT scores and remove any NaN values
    scores = np.asarray(df["HGT_score"].dropna().values, dtype=np.float64)
    logger.info("Analyzing %d HGT scores (after removing NaN)", len(scores))

    # Compute basic statistics
    logger.info("Computing basic statistics...")
    stats = compute_basic_statistics(scores)
    logger.info(
        "Stats: count=%d, mean=%.4f, median=%.4f, q95=%.4f",
        stats.count,
        stats.mean,
        stats.median,
        stats.q95,
    )

    # Fit theoretical distributions
    logger.info("Fitting theoretical distributions...")
    fitted_dists = fit_distributions(scores)

    # Bootstrap analysis
    logger.info("Performing bootstrap analysis for threshold determination...")
    bootstrap_result = bootstrap_quantile(scores, quantile=0.95, n_bootstrap=10000)

    # Create visualizations
    logger.info("Generating visualizations...")
    create_comprehensive_plots(scores, fitted_dists, bootstrap_result, output_dir)

    # Generate report
    logger.info("Generating analysis report...")
    report_file = output_dir / "hgt_score_analysis_report.txt"
    generate_report(stats, fitted_dists, bootstrap_result, report_file)

    # Save thresholds to JSON
    if save_threshold_json:
        threshold_file = output_dir / "hgt_thresholds.json"
        save_thresholds(bootstrap_result, threshold_file)

    # Log final summary
    upper_ci = bootstrap_result.upper_ci
    logger.info("=" * 60)
    logger.info("ANALYSIS COMPLETE - KEY FINDINGS")
    logger.info("=" * 60)
    logger.info("RECOMMENDED HGT THRESHOLD: %.6f", upper_ci)
    logger.info(
        "Pairs with HGT_score > %.6f should be considered as "
        "significantly high and potential candidates for true HGT events.",
        upper_ci,
    )

    # Filter for High Confidence HGT Candidates
    significant_hits = df[df["HGT_score"] > upper_ci]
    output_file = output_dir / "high_confidence_hgt.tsv"
    significant_hits.to_csv(output_file, sep="\t", index=False)
    logger.info(
        "Identified %d high-confidence HGT candidates, saved to %s",
        len(significant_hits),
        output_file,
    )
