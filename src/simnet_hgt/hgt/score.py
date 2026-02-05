"""Core HGT (Horizontal Gene Transfer) score computation.

This module provides the main scoring functions for computing HGT scores
from protein similarity data and phylogenetic distances.
"""

import csv
import json
import logging
import math
import time
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Any

from tqdm import tqdm

from .loaders import (
    get_canonical_genome_key,
    load_genome_to_taxid,
    load_phylogeny_distances,
    load_protein_to_chg,
    load_protein_to_genome,
)
from .models import HgtParams, ScoringMethod, SimilarityStats
from .ranks import build_adjacency, compute_simple_ranks, compute_threshold_ranks
from .scorers import get_scorer

logger = logging.getLogger(__name__)

# Constants
MIN_SAMPLES_FOR_VARIANCE = 2


def compute_similarity_mean(
    hits_file: Path,
    cache_file: Path | None = None,
    has_header: bool = True,
    max_lines: int | None = None,
) -> float:
    """Calculate mean similarity (pident) by streaming the hits file.

    Supports caching to avoid recomputation on subsequent runs.

    Args:
        hits_file: Path to TSV file with protein hits.
            Expected columns: proteinA, proteinB, pident, ...
        cache_file: Optional path to JSON file for caching the result.
            If exists and valid, loads from cache instead of recomputing.
        has_header: Whether the file has a header row. Defaults to True.
        max_lines: Optional limit on lines to process.

    Returns:
        Mean similarity value (pident).
    """
    # Try to load from cache if provided
    if cache_file and cache_file.exists():
        logger.info("Loading similarity mean from cache: %s", cache_file)
        with cache_file.open(encoding="utf-8") as f:
            data = json.load(f)
            return data.get("s_mean", 0.0)

    # Compute from scratch
    total = 0.0
    count = 0

    # Estimate line count for progress bar
    line_count: int | None = None
    try:
        with hits_file.open("rb") as f:
            line_count = sum(1 for _ in f)
    except Exception:
        pass

    show_progress = logger.isEnabledFor(logging.INFO)

    with hits_file.open(encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        first = True

        for i, row in tqdm(
            enumerate(reader),
            total=line_count,
            desc="Computing s_mean",
            disable=not show_progress,
        ):
            if first and has_header:
                first = False
                continue
            if max_lines is not None and i > max_lines:
                break
            if not row or len(row) < 3:  # noqa: PLR2004
                continue
            try:
                pident = float(row[2])
            except ValueError:
                continue
            total += pident
            count += 1

    s_mean = (total / count) if count else 0.0

    # Save to cache if provided
    if cache_file:
        logger.info("Saving similarity mean to cache: %s", cache_file)
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with cache_file.open("w", encoding="utf-8") as f:
            json.dump({"s_mean": s_mean}, f, indent=2)

    return s_mean


def compute_similarity_stats_by_distance(
    hits_data: list[dict],
) -> SimilarityStats:
    """Compute mean and std of similarity for each phylogenetic distance.

    Args:
        hits_data: List of hit dictionaries with 'n_int' and 'similarity' keys.

    Returns:
        Dictionary mapping N_int to (mean, std) tuples.
    """
    # Group similarities by N_int
    by_distance: dict[int, list[float]] = defaultdict(list)
    for hit in hits_data:
        n_int = hit["n_int"]
        if n_int > 0:
            by_distance[n_int].append(hit["similarity"])

    # Compute mean and std for each distance
    stats: SimilarityStats = {}
    for n_int, similarities in by_distance.items():
        if len(similarities) >= MIN_SAMPLES_FOR_VARIANCE:
            mean = sum(similarities) / len(similarities)
            variance = sum((s - mean) ** 2 for s in similarities) / (len(similarities) - 1)
            std = math.sqrt(variance) if variance > 0 else 1.0
            stats[n_int] = (mean, std)
        elif len(similarities) == 1:
            stats[n_int] = (similarities[0], 1.0)  # Use std=1 for single samples

    return stats


def compute_taxon_stats(hits: list[dict[str, Any]]) -> dict[str, tuple[float, float]]:
    """Compute mean and std similarity for each TaxID."""
    tax_values: dict[str, list[float]] = defaultdict(list)
    for h in hits:
        s = h["similarity"]
        tax_values[h["tax_a"]].append(s)
        tax_values[h["tax_b"]].append(s)

    stats = {}
    for tax, vals in tax_values.items():
        if not vals:
            continue
        mean = sum(vals) / len(vals)
        variance = sum((x - mean) ** 2 for x in vals) / len(vals) if len(vals) > 1 else 0.0
        std = math.sqrt(variance)
        # Regularize std slightly to avoid division by zero
        stats[tax] = (mean, max(std, 2.0))
    return stats


def compute_n_int_cdf(hits: dict[tuple, dict]) -> dict[int, float]:
    """Compute cumulative distribution function for N_int values."""
    n_values = sorted([d["n_int"] for d in hits.values()])
    total_n = len(n_values)
    if total_n == 0:
        return {}

    n_cdf: dict[int, float] = {}
    unique_ns = sorted(set(n_values))
    for n in unique_ns:
        rank = bisect_right(n_values, n)
        n_cdf[n] = rank / total_n
    return n_cdf


def _collect_best_hits(  # noqa: PLR0913
    hits_tsv: Path,
    protein_to_chg: dict[str, str],
    protein_to_genome: dict[str, str],
    phylogeny_distances: dict[tuple, int],
    accession_to_taxid_map: dict[str, str],
    scoring_method: ScoringMethod,
) -> dict[tuple, dict]:
    """First pass: collect best hits per TaxID pair within each CHG."""
    # Estimate line count for progress bar
    line_count: int | None = None
    try:
        with hits_tsv.open("rb") as f:
            line_count = sum(1 for _ in f) - 1
    except Exception:
        pass

    best_taxid_hits: dict[tuple, dict] = {}
    show_progress = logger.isEnabledFor(logging.INFO)
    is_ranking = scoring_method in (ScoringMethod.RANKING, ScoringMethod.RANKING_TH)

    with hits_tsv.open(encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)  # Skip header

        for row in tqdm(
            reader,
            total=line_count,
            desc="Computing HGT scores",
            disable=not show_progress,
        ):
            prot_a, prot_b, similarity = row[0], row[1], float(row[2])

            # Map proteins to their respective TaxIDs via their Genomes
            genome_a = protein_to_genome.get(prot_a)
            genome_b = protein_to_genome.get(prot_b)

            if not genome_a or not genome_b:
                continue

            tax_a = accession_to_taxid_map.get(genome_a)
            tax_b = accession_to_taxid_map.get(genome_b)

            if tax_a == tax_b or not tax_a or not tax_b:
                continue

            chg_id = protein_to_chg.get(prot_a)
            if chg_id != protein_to_chg.get(prot_b):
                continue

            # Key for TaxID grouping
            if is_ranking:
                # Directional key for asymmetric scoring
                tax_key = (chg_id, tax_a, tax_b)
            else:
                # Symmetric key for standard scoring
                tax_key = (chg_id, tuple(sorted((tax_a, tax_b))))

            # Keep only the largest similarity value for this pair of taxa
            if (
                tax_key not in best_taxid_hits
                or similarity > best_taxid_hits[tax_key]["similarity"]
            ):
                n_int = phylogeny_distances.get(get_canonical_genome_key(genome_a, genome_b), 0)
                best_taxid_hits[tax_key] = {
                    "prot_a": prot_a,
                    "prot_b": prot_b,
                    "genome_a": genome_a,
                    "genome_b": genome_b,
                    "tax_a": tax_a,
                    "tax_b": tax_b,
                    "n_int": n_int,
                    "similarity": similarity,
                }

    return best_taxid_hits


def _symmetrize_hits(best_taxid_hits: dict[tuple, dict]) -> dict[tuple, dict]:
    """Symmetrize hits for ranking methods to ensure both directions are scored."""
    new_hits = {}
    for key, data in best_taxid_hits.items():
        chg_id, t_a, t_b = key
        reverse_key = (chg_id, t_b, t_a)
        if reverse_key not in best_taxid_hits:
            new_data = data.copy()
            new_data["prot_a"] = data["prot_b"]
            new_data["prot_b"] = data["prot_a"]
            new_data["genome_a"] = data["genome_b"]
            new_data["genome_b"] = data["genome_a"]
            new_data["tax_a"] = data["tax_b"]
            new_data["tax_b"] = data["tax_a"]
            new_hits[reverse_key] = new_data

    result = dict(best_taxid_hits)
    result.update(new_hits)
    logger.info("Added %d reverse hits", len(new_hits))
    return result


def _prepare_scoring_metadata(
    best_taxid_hits: dict[tuple, dict],
    params: HgtParams,
) -> tuple[SimilarityStats, dict, dict, dict]:
    """Prepare all auxiliary data structures needed for scoring."""
    # Compute similarity statistics for z-score methods
    similarity_stats: SimilarityStats = {}
    if params.scoring_method in (
        ScoringMethod.ZSCORE,
        ScoringMethod.MODIFIED,
        ScoringMethod.COMBINED,
        ScoringMethod.EXPONENTIAL,
    ):
        similarity_stats = compute_similarity_stats_by_distance(list(best_taxid_hits.values()))
        logger.info("Computed similarity stats for %d distance levels", len(similarity_stats))

    # Compute taxon stats and N_int CDF for Exponential method
    taxon_stats: dict[str, tuple[float, float]] = {}
    n_cdf: dict[int, float] = {}
    if params.scoring_method == ScoringMethod.EXPONENTIAL:
        taxon_stats = compute_taxon_stats(list(best_taxid_hits.values()))
        logger.info("Computed taxon statistics for %d taxa", len(taxon_stats))
        n_cdf = compute_n_int_cdf(best_taxid_hits)
        logger.info("Computed N_int CDF for adaptive distance weighting")

    # Compute neighbor ranks for Ranking methods
    neighbor_ranks: dict[str, dict[str, dict[str, int]]] = {}
    if params.scoring_method == ScoringMethod.RANKING:
        logger.info("Computing ranks for Ranking method")
        adj = build_adjacency(best_taxid_hits, params.min_phylogenetic_distance)
        neighbor_ranks = compute_simple_ranks(adj)
    elif params.scoring_method == ScoringMethod.RANKING_TH:
        logger.info("Computing threshold ranks for Ranking TH method")
        adj = build_adjacency(best_taxid_hits, params.min_phylogenetic_distance)
        neighbor_ranks = compute_threshold_ranks(adj, params.rank_threshold)

    return similarity_stats, taxon_stats, n_cdf, neighbor_ranks


def compute_hgt_scores(
    hits_tsv: Path,
    file_paths: dict[str, Path],
    params: HgtParams,
) -> Path:
    """Compute HGT scores from protein hits data.

    Processes protein similarity hits, maps them to taxa, and computes
    HGT scores. Filters for the largest similarity value per TaxID pair
    within each CHG as required by the methodology.

    Args:
        hits_tsv: Path to TSV file with protein hits.
            Expected columns: proteinA, proteinB, pident, ...
        file_paths: Dictionary with paths to required data files:
            - "protein_to_genome": CSV mapping proteins to genomes
            - "protein_to_chg": TSV mapping proteins to CHGs
            - "phylogeny_distances": TSV with pairwise phylogenetic distances
            - "genome_to_taxid": CSV mapping genomes to taxids
            - "output": Directory for output files
        params: HGT computation parameters.

    Returns:
        Tuple of (Path to output file, computation time in seconds)
    """
    # Load required mappings
    protein_to_chg = load_protein_to_chg(file_paths["protein_to_chg"])
    protein_to_genome = load_protein_to_genome(file_paths["protein_to_genome"])
    phylogeny_distances = load_phylogeny_distances(file_paths["phylogeny_distances"])
    accession_to_taxid_map = load_genome_to_taxid(file_paths["genome_to_taxid"])

    # Create output directory and file path
    output_dir = file_paths["output"]
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "hgt_scores.tsv"

    start_time = time.perf_counter()

    # Collect best hits
    best_taxid_hits = _collect_best_hits(
        hits_tsv,
        protein_to_chg,
        protein_to_genome,
        phylogeny_distances,
        accession_to_taxid_map,
        params.scoring_method,
    )

    # Symmetrize hits for Ranking methods
    # Symmetrize hits for Ranking methods
    is_ranking = params.scoring_method in (ScoringMethod.RANKING, ScoringMethod.RANKING_TH)
    if is_ranking:
        logger.info("Symmetrizing hits for Ranking method")
        best_taxid_hits = _symmetrize_hits(best_taxid_hits)

    # Compute detailed statistics and ranks
    (
        similarity_stats,
        taxon_stats,
        n_cdf,
        neighbor_ranks,
    ) = _prepare_scoring_metadata(best_taxid_hits, params)

    # Get scorer function
    scorer = get_scorer(params.scoring_method)

    # Find max N_int for combined score normalization
    filtered_hits = [
        d for d in best_taxid_hits.values() if d["n_int"] >= params.min_phylogenetic_distance
    ]
    n_int_max = max((d["n_int"] for d in filtered_hits), default=10)

    # Calculate scores and build results
    results = []
    for tax_key, data in best_taxid_hits.items():
        if data["n_int"] < params.min_phylogenetic_distance:
            continue

        # Default rank for non-ranking methods
        rank = 0
        if is_ranking:
            chg_id = tax_key[0]
            rank = neighbor_ranks.get(chg_id, {}).get(data["tax_a"], {}).get(data["tax_b"], 0)

        # Build kwargs for scorer
        scorer_kwargs = {
            "n_int": data["n_int"],
            "similarity": data["similarity"],
            "n_int_mean": params.n_int_mean,
            "s_mean": params.hits_mean,
            "similarity_stats": similarity_stats,
            "taxon_stats": taxon_stats,
            "tax_a": data["tax_a"],
            "tax_b": data["tax_b"],
            "n_cdf": n_cdf,
            "n_int_max": n_int_max,
            "rank": rank,
        }

        score = scorer(**scorer_kwargs)

        results.append(
            [
                data["prot_a"],
                data["prot_b"],
                tax_key[0],
                data["tax_a"],
                data["tax_b"],
                data["n_int"],
                f"{data['similarity']:.1f}",
                rank,
                score,
            ]
        )

    # Sort results by HGT_score descending and write output
    # Sort results by HGT_score descending and write output
    results.sort(key=lambda x: x[8], reverse=True)

    end_time = time.perf_counter()

    _write_hgt_results(output_file, results)

    logger.info("Wrote HGT scores to %s (method: %s)", output_file, params.scoring_method.value)
    return output_file, end_time - start_time


def _write_hgt_results(output_file: Path, results: list) -> None:
    """Write HGT scores to TSV file."""
    with output_file.open("w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(
            [
                "protein_A",
                "protein_B",
                "CHG",
                "TaxID_A",
                "TaxID_B",
                "N_int",
                "s_AB",
                "Rank",
                "HGT_score",
            ]
        )
        writer.writerows(results)
