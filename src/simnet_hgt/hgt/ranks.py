"""Rank computation for HGT scoring.

This module provides functions to compute neighbor ranks for ranking-based
HGT scoring methods.
"""

import logging
from collections import defaultdict
from typing import Any

logger = logging.getLogger(__name__)

# Type alias: chg_id -> tax_id -> neighbor_tax_id -> rank
NeighborRanks = dict[str, dict[str, dict[str, int]]]


def build_adjacency(
    hits: dict[tuple, dict[str, Any]],
    min_phylogenetic_distance: int,
) -> dict[str, dict[str, dict[str, float]]]:
    """Build adjacency map from hits data.

    Creates a mapping: chg_id -> tax_id -> neighbor_tax_id -> max_similarity

    Args:
        hits: Dictionary of (chg_id, tax_a, tax_b) -> hit data.
        min_phylogenetic_distance: Minimum N_int to include a hit.

    Returns:
        Nested dictionary for adjacency lookups.
    """
    adj: dict[str, dict[str, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))

    for key, data in hits.items():
        if data["n_int"] < min_phylogenetic_distance:
            continue

        chg_id = key[0]
        t_a, t_b = data["tax_a"], data["tax_b"]
        sim = data["similarity"]

        # Update max similarity for t_a -> t_b
        if t_b not in adj[chg_id][t_a] or sim > adj[chg_id][t_a][t_b]:
            adj[chg_id][t_a][t_b] = sim

        # Update max similarity for t_b -> t_a (symmetric)
        if t_a not in adj[chg_id][t_b] or sim > adj[chg_id][t_b][t_a]:
            adj[chg_id][t_b][t_a] = sim

    return adj


def compute_simple_ranks(
    adj: dict[str, dict[str, dict[str, float]]],
) -> NeighborRanks:
    """Compute sequential ranks based on similarity ordering.

    Each neighbor is assigned a rank based on its position when sorted
    by similarity (descending). Rank 1 = highest similarity.

    Args:
        adj: Adjacency map from build_adjacency().

    Returns:
        Nested dictionary mapping chg_id -> tax_id -> neighbor_tax_id -> rank.
    """
    neighbor_ranks: NeighborRanks = defaultdict(lambda: defaultdict(dict))

    for chg_id, tax_map in adj.items():
        for tax_id, neighbor_map in tax_map.items():
            # Sort by similarity descending
            neighbors = sorted(neighbor_map.items(), key=lambda x: x[1], reverse=True)
            # Assign sequential ranks
            for rank_idx, (neighbor_tax, _) in enumerate(neighbors):
                neighbor_ranks[chg_id][tax_id][neighbor_tax] = rank_idx + 1

    return neighbor_ranks


def compute_threshold_ranks(
    adj: dict[str, dict[str, dict[str, float]]],
    threshold: float,
) -> NeighborRanks:
    """Compute grouped ranks based on similarity threshold.

    Neighbors with similarity differences <= threshold share the same rank.
    Only increments rank when the gap between consecutive similarities
    exceeds the threshold.

    Args:
        adj: Adjacency map from build_adjacency().
        threshold: Maximum similarity difference for same-rank grouping.

    Returns:
        Nested dictionary mapping chg_id -> tax_id -> neighbor_tax_id -> rank.
    """
    neighbor_ranks: NeighborRanks = defaultdict(lambda: defaultdict(dict))

    for chg_id, tax_map in adj.items():
        for tax_id, neighbor_map in tax_map.items():
            if not neighbor_map:
                continue

            # Sort by similarity descending
            neighbors = sorted(neighbor_map.items(), key=lambda x: x[1], reverse=True)

            # First neighbor always gets rank 1
            current_rank = 1
            neighbor_ranks[chg_id][tax_id][neighbors[0][0]] = current_rank

            # Iterate from second element
            for i in range(1, len(neighbors)):
                prev_sim = neighbors[i - 1][1]
                curr_tax, curr_sim = neighbors[i]

                # If diff > threshold, increment rank
                if (prev_sim - curr_sim) > threshold:
                    current_rank += 1

                neighbor_ranks[chg_id][tax_id][curr_tax] = current_rank

    return neighbor_ranks
