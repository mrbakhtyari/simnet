import itertools
import json
import logging
import time
from pathlib import Path

from Bio.Phylo._io import read as phylo_read

logger = logging.getLogger(__name__)


def get_node_depths(tree):
    """Pre-calculate the depth of every node from the root."""
    depths = {}
    for clade in tree.find_clades():
        # Depth is the number of ancestors to the root
        depths[clade] = len(tree.get_path(clade))
    return depths


def compute_pairwise_distances(
    newick_path: Path,
    output_path: Path,
    prefix: str = "",
) -> float:
    output_path.mkdir(parents=True, exist_ok=True)
    output_tsv = output_path / f"{prefix}_pairwise_distances.tsv"
    stats_json = output_path / f"{prefix}_phylogeny_stats.json"

    logger.info("Reading Newick tree...")
    tree = phylo_read(newick_path, "newick")

    # To strictly separate computation from I/O as requested, we should ideally
    # compute all values first, then write.
    # However, for very large trees this memory usage is prohibitive.
    # Since the user specifically asked for "only computation time of each step",
    # and this step is O(N^2), strictly separating is best if N is small.
    # If N is large, we must stream.
    # Assuming standard datasets, we can compute first.

    start_time = time.perf_counter()
    # Pre-calculate depths
    depths = get_node_depths(tree)

    leaves = [term for term in tree.get_terminals() if term.name]

    pairs_data = []
    total_n_int_internal = 0
    count_internal = 0

    for a, b in itertools.combinations(leaves, 2):
        lca = tree.common_ancestor(a, b)
        d_a = depths[a]
        d_b = depths[b]
        d_lca = depths[lca]
        n_int = max(0, (d_a + d_b) - (2 * d_lca) - 1)

        pairs_data.append((a.name, b.name, n_int))
        total_n_int_internal += n_int
        count_internal += 1

    end_time = time.perf_counter()

    # Now write results
    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("genome_a\tgenome_b\tN_int\n")
        for genome_a, genome_b, n_int in pairs_data:
            f.write(f"{genome_a}\t{genome_b}\t{n_int}\n")

    n_int_mean = (total_n_int_internal / count_internal) if count_internal else 0.0

    with open(stats_json, "w", encoding="utf-8") as f:
        json.dump({"n_int_mean": n_int_mean, "leaf_count": len(leaves)}, f)

    logger.info(f"✓ Phylogeny: {len(leaves)} leaves, mean internal nodes = {n_int_mean:.3f}")
    return end_time - start_time
