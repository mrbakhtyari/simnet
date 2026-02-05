import itertools
import logging
from collections import Counter, defaultdict
from pathlib import Path

from .io import (
    get_category_color,
    iter_cluster_members,
    load_protein_metadata,
    write_edges,
    write_nodes,
)
from .network import NetworkEdge, NetworkNode

logger = logging.getLogger(__name__)


def _build_chg_to_genomes(
    chg_membership_tsv: Path,
    protein_to_genome: dict[str, str],
) -> dict[str, set[str]]:
    """Build mapping from CHG to set of genomes containing it.

    Args:
        chg_membership_tsv: CHG membership file
        protein_to_genome: Protein -> genome mapping

    Returns:
        Dictionary mapping CHG ID to set of genome IDs
    """
    chg_to_genomes: dict[str, set[str]] = defaultdict(set)
    for protein_id, chg_id in iter_cluster_members(chg_membership_tsv):
        genome = protein_to_genome.get(protein_id)
        if genome is not None:
            chg_to_genomes[chg_id].add(genome)
    return chg_to_genomes


def build_genome_network(
    chg_membership_tsv: Path,
    protein_metadata_csv: Path,
    output_dir: Path,
) -> None:
    p2g, g2t = load_protein_metadata(protein_metadata_csv)
    chg_to_genomes = _build_chg_to_genomes(chg_membership_tsv, p2g)
    nodes, edges = _build_genome_similarity_from_chg_map(chg_to_genomes, g2t)
    nodes_csv = output_dir / "genome_nodes.csv"
    edges_csv = output_dir / "genome_edges.csv"

    write_nodes(nodes_csv, nodes, include_reduction_fields=False)
    write_edges(edges_csv, edges)

    logger.info(
        f"Genome similarity network built with {len(nodes)} nodes and {len(edges)} edges. "
        f"Nodes written to {nodes_csv}, edges to {edges_csv}."
    )


def _build_genome_similarity_from_chg_map(
    chg_to_genomes: dict[str, set[str]],
    genome_to_category: dict[str, str],
) -> tuple[list[NetworkNode], list[NetworkEdge]]:
    """Build genome similarity network nodes and edges from CHG-to-genomes mapping.

    Args:
        chg_to_genomes: Mapping from CHG ID to set of genome IDs
        genome_to_category: Genome -> category mapping

    Returns:
        Tuple of (nodes, edges) lists
    """
    # Weighted edges: for each CHG, add 1 to each genome pair
    edge_weights: Counter[tuple[str, str]] = Counter()
    genomes_seen: set[str] = set()
    for genomes in chg_to_genomes.values():
        g_list = sorted(genomes)
        genomes_seen.update(g_list)
        for a, b in itertools.combinations(g_list, 2):
            edge = (a, b) if a < b else (b, a)
            edge_weights[edge] += 1

    # Create node objects
    nodes = [
        NetworkNode(
            node_id=g,
            node_type="genome",
            category=genome_to_category.get(g, "unknown"),
            color=get_category_color(genome_to_category.get(g, "unknown")),
        )
        for g in sorted(genomes_seen)
    ]

    # Create edge objects
    edges = [
        NetworkEdge(source=a, target=b, weight=wt) for (a, b), wt in sorted(edge_weights.items())
    ]

    return nodes, edges
