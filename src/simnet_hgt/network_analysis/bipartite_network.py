"""Functions for building biological networks from CHG and genome data.

This module provides functions to create two types of networks:
1. Genome-Genome network weighted by shared CHGs
2. Genome-CHG bipartite network with protein counts as weights
"""

import logging
import time
from collections import Counter
from pathlib import Path

from .io import (
    CATEGORY_COLORS,
    iter_cluster_members,
    load_protein_metadata,
    write_edges,
    write_nodes,
)
from .network import NetworkEdge, NetworkNode

logger = logging.getLogger(__name__)


def build_bipartite_network(
    connected_components_tsv: Path,
    protein_metadata_csv: Path,
    output_dir: Path,
) -> float:
    """Build biological networks from CHG and protein metadata.

    Args:
        connected_components_tsv: CHG membership file
        protein_metadata_csv: Unified protein metadata file (protein_id, genome_accession, taxon)
        output_dir: Output directory for nodes and edges CSV files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load mappings from unified metadata file
    p2g, g2t = load_protein_metadata(protein_metadata_csv)

    # Build genome similarity network (if requested)
    # Always build bipartite network
    # Load all inputs first

    # For strict computation timing, we should measure the logic inside _compute_genome_chg_counts
    # and _build_bipartite_from_counts.
    # However, _compute_genome_chg_counts iterates the file.
    # To truly isolate, we'd need to load the file into memory first.
    # Assuming iter_cluster_members is efficient streaming, we will load memberships into memory.

    # 1. Load data
    memberships = list(iter_cluster_members(connected_components_tsv))

    start_time = time.perf_counter()
    # Logic: process memberships -> counts
    counts: Counter[tuple[str, str]] = Counter()
    chgs_seen: set[str] = set()
    genomes_seen: set[str] = set()

    for protein_id, chg_id in memberships:
        genome = p2g.get(protein_id)
        if genome is None:
            continue
        counts[(genome, chg_id)] += 1
        chgs_seen.add(chg_id)
        genomes_seen.add(genome)

    genome_chg_counts = (counts, chgs_seen, genomes_seen)

    nodes, edges = _build_bipartite_from_counts(genome_chg_counts, g2t)
    end_time = time.perf_counter()

    nodes_csv = output_dir / "bipartite_network_nodes.csv"
    edges_csv = output_dir / "bipartite_network_edges.csv"

    write_nodes(nodes_csv, nodes, include_reduction_fields=False)
    write_edges(edges_csv, edges)

    logger.info(
        f"Bipartite network built with {len(nodes)} nodes and {len(edges)} edges. "
        f"Nodes written to {nodes_csv}, edges to {edges_csv}."
    )
    return end_time - start_time


def _compute_genome_chg_counts(
    connected_components_tsv: Path,
    protein_to_genome: dict[str, str],
) -> tuple[Counter[tuple[str, str]], set[str], set[str]]:
    """Compute protein counts per (genome, CHG) pair.

    Args:
        connected_components_tsv: CHG membership file
        protein_to_genome: Protein -> genome mapping

    Returns:
        Tuple of (counts, chgs_seen, genomes_seen)
    """
    counts: Counter[tuple[str, str]] = Counter()
    chgs_seen: set[str] = set()
    genomes_seen: set[str] = set()

    for protein_id, chg_id in iter_cluster_members(connected_components_tsv):
        genome = protein_to_genome.get(protein_id)
        if genome is None:
            continue
        counts[(genome, chg_id)] += 1
        chgs_seen.add(chg_id)
        genomes_seen.add(genome)

    return counts, chgs_seen, genomes_seen


def _build_bipartite_from_counts(
    genome_chg_counts: tuple[Counter[tuple[str, str]], set[str], set[str]],
    genome_to_category: dict[str, str],
) -> tuple[list[NetworkNode], list[NetworkEdge]]:
    """Build bipartite network nodes and edges from pre-computed counts.

    Filters out CHGs that connect to only one genome (trivial clusters).

    Args:
        genome_chg_counts: Tuple of (counts, chgs_seen, genomes_seen) from _compute_genome_chg_counts
        genome_to_category: Genome -> category mapping

    Returns:
        Tuple of (nodes, edges) lists
    """
    counts, chgs_seen, _genomes_seen = genome_chg_counts

    # Count how many genomes each CHG connects to
    chg_genome_count: Counter[str] = Counter()
    for _genome, chg in counts:
        chg_genome_count[chg] += 1

    # Filter: keep only CHGs with 2+ genome connections (non-trivial clusters)
    min_genomes = 2
    valid_chgs = {chg for chg, count in chg_genome_count.items() if count >= min_genomes}
    filtered_chgs_count = len(chgs_seen) - len(valid_chgs)

    # Filter counts to only include valid CHGs
    filtered_counts = {(g, chg): wt for (g, chg), wt in counts.items() if chg in valid_chgs}

    # Update genomes to only include those that have connections to valid CHGs
    filtered_genomes = {g for g, _chg in filtered_counts}

    # Create nodes: genomes + CHGs
    nodes: list[NetworkNode] = []

    # Add genome nodes (only those with valid CHG connections)
    for g in sorted(filtered_genomes):
        cat = genome_to_category.get(g, "unknown")
        nodes.append(
            NetworkNode(
                node_id=g,
                node_type="genome",
                category=cat,
                color=CATEGORY_COLORS.get(cat, "#7f7f7f"),
            )
        )

    # Add CHG nodes (only valid ones)
    for chg in sorted(valid_chgs):
        nodes.append(NetworkNode(node_id=chg, node_type="CHG", category="", color="#bbbbbb"))

    # Create edges: genome <-> CHG with weight (only for valid CHGs)
    edges = [
        NetworkEdge(source=g, target=chg, weight=wt) for (g, chg), wt in filtered_counts.items()
    ]

    logger.info(
        f"Filtered out {filtered_chgs_count} CHGs with only 1 genome connection. "
        f"Kept {len(valid_chgs)} CHGs with 2+ genome connections."
    )

    return nodes, edges
