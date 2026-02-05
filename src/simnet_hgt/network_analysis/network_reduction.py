"""Functions for reducing bipartite graphs by merging twin CHG nodes.

Twin nodes are CHG nodes with identical genome neighborhoods. This module:
1. Detects twin CHG nodes
2. Merges them into representative nodes
3. Identifies articulation points (candidate genetic public goods)
"""

import csv
import logging
import time
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx

from .io import (
    iter_edges_simple,
    read_nodes,
    write_edges,
    write_nodes,
)
from .network import NetworkEdge, NetworkNode

logger = logging.getLogger(__name__)


def reduce_bipartite_network(
    bipartite_nodes_csv: Path,
    bipartite_edges_csv: Path,
    connected_components_tsv: Path,
) -> list[str]:
    """Reduce bipartite graph by merging twin CHG nodes.

    Args:
        bipartite_nodes_csv: Input nodes
        bipartite_edges_csv: Input edges
        chg_membership_tsv: Original CHG membership

    Returns:
        Tuple of (list of reduced node IDs, computation time in seconds)
    """

    # Derive output paths from input directory
    input_dir = bipartite_nodes_csv.parent
    reduced_nodes_csv = input_dir / "reduced_bipartite_nodes.csv"
    reduced_edges_csv = input_dir / "reduced_bipartite_edges.csv"
    twin_groups_tsv = input_dir / "chg_twin_groups.tsv"
    reduced_chg_membership_tsv = input_dir / "reduced_chg_membership.tsv"

    # Read node metadata
    nodes = read_nodes(bipartite_nodes_csv)
    node_meta = {n.node_id: n for n in nodes}

    # Build adjacency for CHGs: chg -> Counter(genome -> weight)
    chg_adj: dict[str, Counter[str]] = defaultdict(Counter)
    genome_nodes: set[str] = set()
    chg_nodes: set[str] = set()

    for a, b, w in iter_edges_simple(bipartite_edges_csv):
        a_meta = node_meta.get(a)
        b_meta = node_meta.get(b)
        if not a_meta or not b_meta:
            continue

        # Determine which is genome and which is CHG
        if a_meta.node_type == "genome" and b_meta.node_type == "CHG":
            genome, chg = a, b
        elif a_meta.node_type == "CHG" and b_meta.node_type == "genome":
            genome, chg = b, a
        else:
            continue

        chg_nodes.add(chg)
        genome_nodes.add(genome)
        chg_adj[chg][genome] += w

    # Start timing pure reduction logic
    start_time = time.perf_counter()

    # Compute twin groups and representatives
    chg_to_rep, rep_to_members = _identify_twin_chgs(chg_adj)

    # Build reduced edges: aggregate weights per (genome, rep)
    agg = _aggregate_edges(chg_adj, chg_to_rep)

    # Build NetworkX graph to compute articulation points
    G = nx.Graph()
    for (g, rep), w in agg.items():
        G.add_edge(g, rep, weight=w)

    # Mark articulation points (public goods candidates)
    artics = set(nx.articulation_points(G))

    # Create reduced nodes
    reduced_nodes: list[NetworkNode] = []

    # Add genome nodes
    for g in sorted(genome_nodes):
        m = node_meta.get(g)
        if m:
            reduced_nodes.append(
                NetworkNode(
                    node_id=g,
                    node_type="genome",
                    category=m.category,
                    color=m.color,
                    twin_group="",
                    is_articulation=g in artics,
                )
            )

    # Add CHG representative nodes
    for rep, _members in rep_to_members.items():
        # If rep is an original CHG (singleton), keep original color/meta
        if rep in node_meta:
            m = node_meta[rep]
            color = m.color
        else:
            color = "#bbbbbb"

        group_label = rep if rep.startswith("CHG_TWIN_") else ""
        reduced_nodes.append(
            NetworkNode(
                node_id=rep,
                node_type="CHG",
                category="",
                color=color,
                twin_group=group_label,
                is_articulation=rep in artics,
            )
        )

    # Create reduced edges
    reduced_edges = [
        NetworkEdge(source=g, target=rep, weight=wt) for (g, rep), wt in sorted(agg.items())
    ]

    end_time = time.perf_counter()

    # Write reduced nodes and edges
    write_nodes(reduced_nodes_csv, reduced_nodes, include_reduction_fields=True)
    write_edges(reduced_edges_csv, reduced_edges)

    # Write twin groups for traceability
    _write_twin_groups(twin_groups_tsv, rep_to_members)

    # Write reduced CHG membership
    _write_reduced_chg_membership(
        connected_components_tsv,
        reduced_chg_membership_tsv,
        chg_to_rep,
    )
    logger.info(
        f"Graph reduction complete. Reduced nodes written to {reduced_nodes_csv}, "
        f"edges to {reduced_edges_csv}, twin groups to {twin_groups_tsv}, "
        f"and reduced CHG membership to {reduced_chg_membership_tsv}."
    )

    # Return list of genome node IDs only (excluding CHG nodes)
    return [
        node.node_id for node in reduced_nodes if node.node_type == "genome"
    ], end_time - start_time


def _identify_twin_chgs(
    chg_adj: dict[str, Counter[str]],
) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Compute twin groups based on identical genome neighborhoods.

    Args:
        chg_adj: CHG -> genome neighbors adjacency

    Returns:
        (chg_to_rep, rep_to_members) mappings
    """
    groups: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for chg, counter in chg_adj.items():
        key = tuple(sorted(counter.keys()))
        groups[key].append(chg)

    chg_to_rep: dict[str, str] = {}
    rep_to_members: dict[str, list[str]] = {}
    group_idx = 0

    for members in groups.values():
        if len(members) == 1:
            # Singleton: CHG is its own representative
            chg_to_rep[members[0]] = members[0]
            rep_to_members[members[0]] = [members[0]]
        else:
            # Twin group: create a representative
            group_idx += 1
            rep = f"CHG_TWIN_{group_idx}"
            rep_to_members[rep] = list(members)
            for m in members:
                chg_to_rep[m] = rep

    return chg_to_rep, rep_to_members


def _aggregate_edges(
    chg_adj: dict[str, Counter[str]], chg_to_rep: dict[str, str]
) -> dict[tuple[str, str], int]:
    """Aggregate edges by representative CHG nodes.

    Args:
        chg_adj: Original CHG adjacency
        chg_to_rep: CHG to representative mapping

    Returns:
        Aggregated edges as (genome, rep) -> weight
    """
    agg: dict[tuple[str, str], int] = defaultdict(int)
    for chg, neigh in chg_adj.items():
        rep = chg_to_rep[chg]
        for g, w in neigh.items():
            agg[(g, rep)] += w
    return dict(agg)


def _write_twin_groups(twin_groups_tsv: Path, rep_to_members: dict[str, list[str]]) -> None:
    """Write twin group membership to TSV."""
    twin_groups_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(twin_groups_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["twin_group", "member_chg"])
        for rep, members in sorted(rep_to_members.items()):
            if rep.startswith("CHG_TWIN_"):
                for m in members:
                    w.writerow([rep, m])


def _write_reduced_chg_membership(
    chg_membership_path: Path,
    reduced_membership_path: Path,
    chg_to_rep: dict[str, str],
) -> None:
    """Map proteins to their reduced CHG representative."""
    reduced_membership_path.parent.mkdir(parents=True, exist_ok=True)

    with open(chg_membership_path, encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        rows = list(reader)

    with open(reduced_membership_path, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["protein_id", "chg_id", "original_chg_id"])

        for row in rows:
            protein_id = row["protein_id"]
            original_chg = row["chg_id"]
            reduced_chg = chg_to_rep.get(original_chg, original_chg)
            writer.writerow([protein_id, reduced_chg, original_chg])
