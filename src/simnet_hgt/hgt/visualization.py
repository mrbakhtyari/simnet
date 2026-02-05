"""Visualization functions for HGT score analysis.

This module provides plotting functions for visualizing HGT score
distributions, statistical analysis results, and threshold determination.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.axes import Axes
from numpy.typing import NDArray

from .analysis import get_best_fit
from .models import BootstrapResult, DistributionFit

logger = logging.getLogger(__name__)

# Set style for better visualizations
sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (15, 10)
plt.rcParams["font.size"] = 10


def plot_empirical(
    ax: Axes,
    scores: NDArray[np.floating],
    point_estimate: float,
    upper_ci: float,
    q99: float,
) -> None:
    """Plot empirical distribution histogram with threshold lines.

    Args:
        ax: Matplotlib axes to plot on.
        scores: Array of HGT score values.
        point_estimate: 95th percentile point estimate.
        upper_ci: Upper confidence interval bound (threshold).
        q99: 99th percentile value.
    """
    ax.hist(scores, bins=100, alpha=0.7, color="skyblue", edgecolor="black")
    ax.axvline(
        upper_ci, color="red", linestyle="--", linewidth=2.5, label=f"Threshold: {upper_ci:.4f}"
    )
    ax.axvline(
        q99, color="orange", linestyle="--", linewidth=2.5, label=f"99th percentile: {q99:.4f}"
    )
    ax.axvline(
        point_estimate,
        color="green",
        linestyle=":",
        linewidth=2,
        label=f"95th percentile: {point_estimate:.4f}",
    )
    ax.set_xlabel("HGT_score", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.legend(fontsize=12, loc="upper right")
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.3)


def plot_theoretical(
    ax: Axes,
    scores: NDArray[np.floating],
    fitted_dists: list[DistributionFit],
) -> None:
    """Plot histogram with best theoretical distribution overlay.

    Args:
        ax: Matplotlib axes to plot on.
        scores: Array of HGT score values.
        fitted_dists: List of fitted distribution results.
    """
    best_fit = get_best_fit(fitted_dists)

    _n_hist, bins_hist, _patches = ax.hist(
        scores, bins=100, alpha=0.6, color="lightgray", edgecolor="black", label="Empirical"
    )

    # Scale the theoretical PDF to match histogram counts
    x_range = np.linspace(float(scores.min()), float(scores.max()), 1000)
    bin_width = bins_hist[1] - bins_hist[0]
    pdf_scaled = best_fit.distribution.pdf(x_range, *best_fit.params) * len(scores) * bin_width

    ax.plot(
        x_range,
        pdf_scaled,
        "r-",
        linewidth=2.5,
        label=f"{best_fit.name.capitalize()} (KS={best_fit.ks_statistic:.4f})",
    )

    ax.set_xlabel("HGT_score", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(
        f"Best Fit: {best_fit.name.capitalize()} Distribution", fontsize=14, fontweight="bold"
    )
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)


def plot_cdf(
    ax: Axes,
    scores: NDArray[np.floating],
    fitted_dists: list[DistributionFit],
    point_estimate: float,
) -> None:
    """Plot Cumulative Distribution Function comparison.

    Args:
        ax: Matplotlib axes to plot on.
        scores: Array of HGT score values.
        fitted_dists: List of fitted distribution results.
        point_estimate: 95th percentile point estimate.
    """
    best_fit = get_best_fit(fitted_dists)

    sorted_scores = np.sort(scores)
    cdf = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    ax.plot(sorted_scores, cdf, linewidth=2.5, label="Empirical CDF", color="navy")

    # Add best fit theoretical CDF
    x_range = np.linspace(float(scores.min()), float(scores.max()), 1000)
    ax.plot(
        x_range,
        best_fit.distribution.cdf(x_range, *best_fit.params),
        "r--",
        linewidth=2.5,
        label=f"{best_fit.name.capitalize()} CDF",
    )

    ax.axhline(0.95, color="orange", linestyle=":", linewidth=1.5, alpha=0.7)
    ax.axvline(point_estimate, color="orange", linestyle=":", linewidth=1.5, alpha=0.7)
    ax.set_xlabel("HGT_score", fontsize=12)
    ax.set_ylabel("Cumulative Probability", fontsize=12)
    ax.set_title("Cumulative Distribution Function", fontsize=14, fontweight="bold")
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)


def plot_box(
    ax: Axes,
    scores: NDArray[np.floating],
    point_estimate: float,
    upper_ci: float,
    q99: float,
) -> None:
    """Plot box plot with threshold lines.

    Args:
        ax: Matplotlib axes to plot on.
        scores: Array of HGT score values.
        point_estimate: 95th percentile point estimate.
        upper_ci: Upper confidence interval bound.
        q99: 99th percentile value.
    """
    ax.boxplot(
        scores, vert=True, patch_artist=True, boxprops=dict(facecolor="lightblue", alpha=0.7)
    )
    ax.axhline(
        point_estimate,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"95th percentile: {point_estimate:.4f}",
    )
    ax.axhline(
        upper_ci, color="blue", linestyle="--", linewidth=2, label=f"Upper CI: {upper_ci:.4f}"
    )
    ax.axhline(
        q99, color="orange", linestyle="--", linewidth=2, label=f"99th percentile: {q99:.4f}"
    )
    ax.set_ylabel("HGT_score", fontsize=14)
    ax.legend(fontsize=10, loc="center right", bbox_to_anchor=(1.0, 0.35))
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.3)


def plot_violin(
    ax: Axes,
    scores: NDArray[np.floating],
    point_estimate: float,
    upper_ci: float,
) -> None:
    """Plot violin plot showing score distribution.

    Args:
        ax: Matplotlib axes to plot on.
        scores: Array of HGT score values.
        point_estimate: 95th percentile point estimate.
        upper_ci: Upper confidence interval bound.
    """
    ax.violinplot([scores], positions=[1], showmeans=True, showmedians=True)
    ax.axhline(point_estimate, color="red", linestyle="--", linewidth=2)
    ax.axhline(upper_ci, color="blue", linestyle="--", linewidth=2)
    ax.set_ylabel("HGT_score", fontsize=14)
    ax.set_xticks([1])
    ax.set_xticklabels(["HGT_scores"], fontsize=14)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, alpha=0.3)


def plot_bootstrap(
    ax: Axes,
    scores: NDArray[np.floating],
    upper_ci: float,
) -> None:
    """Plot bootstrap distribution across different percentiles.

    Args:
        ax: Matplotlib axes to plot on.
        scores: Array of HGT score values.
        upper_ci: Upper confidence interval bound (threshold).
    """
    logger.debug("Calculating bootstrap distributions for multiple percentiles...")
    percentiles = [90, 95, 99, 99.5]
    bootstrap_results_mult: list[list[float]] = []

    for p in percentiles:
        boot_quantiles: list[float] = []
        for _ in range(1000):  # Quick bootstrap for comparison
            bootstrap_sample = np.random.choice(scores, size=len(scores), replace=True)
            boot_quantiles.append(float(np.percentile(bootstrap_sample, p)))
        bootstrap_results_mult.append(boot_quantiles)

    positions = list(range(len(percentiles)))
    bp = ax.boxplot(bootstrap_results_mult, positions=positions, patch_artist=True)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{p}%" for p in percentiles])

    for patch in bp["boxes"]:
        patch.set_facecolor("lightcoral")
        patch.set_alpha(0.7)

    # Add actual percentiles
    actual_percentiles = [float(np.percentile(scores, p)) for p in percentiles]
    ax.scatter(
        positions,
        actual_percentiles,
        color="red",
        s=100,
        zorder=3,
        label="Actual percentile",
        marker="D",
    )

    ax.axhline(
        upper_ci,
        color="green",
        linestyle="--",
        linewidth=2,
        alpha=0.7,
        label=f"Threshold: {upper_ci:.4f}",
    )

    ax.set_xlabel("Percentile", fontsize=14)
    ax.set_ylabel("HGT_score", fontsize=14)
    ax.legend(fontsize=12, loc="upper left")
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.3)


def create_comprehensive_plots(
    scores: NDArray[np.floating],
    fitted_dists: list[DistributionFit],
    bootstrap_result: BootstrapResult,
    output_dir: Path,
) -> Path:
    """Create comprehensive visualization plots and save to file.

    Generates a 2x2 grid of plots including:
    - (a) Empirical distribution with thresholds
    - (b) Box plot with thresholds
    - (c) Violin plot
    - (d) Bootstrap distribution across percentiles

    Args:
        scores: Array of HGT score values.
        fitted_dists: List of fitted distribution results (unused but kept for API).
        bootstrap_result: Bootstrap analysis results.
        output_dir: Directory to save the output plot.

    Returns:
        Path to the saved plot file.
    """
    point_est = bootstrap_result.point_estimate
    upper_ci = bootstrap_result.upper_ci
    q99 = float(np.percentile(scores, 99))

    _fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) Empirical distribution
    plot_empirical(axes[0, 0], scores, point_est, upper_ci, q99)
    axes[0, 0].set_title("(a)", fontsize=20, fontweight="bold", loc="left")

    # (b) Box plot
    plot_box(axes[0, 1], scores, point_est, upper_ci, q99)
    axes[0, 1].set_title("(b)", fontsize=20, fontweight="bold", loc="left")

    # (c) Violin plot
    plot_violin(axes[1, 0], scores, point_est, upper_ci)
    axes[1, 0].set_title("(c)", fontsize=20, fontweight="bold", loc="left")

    # (d) Bootstrap distribution
    plot_bootstrap(axes[1, 1], scores, upper_ci)
    axes[1, 1].set_title("(d)", fontsize=20, fontweight="bold", loc="left")

    plt.tight_layout()

    output_path = output_dir / "hgt_score_analysis.pdf"
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    logger.info("Saved comprehensive plot to %s", output_path)
    plt.close()

    return output_path
