import logging
from pathlib import Path

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)


def detect_communities(hits_tsv: Path, chg_membership_tsv: Path, min_bitscore: int = 100) -> None:
    """
    Detect protein communities using Louvain community detection algorithm.

    This function builds a protein similarity network from alignment hits and
    identifies communities (clusters of highly similar proteins) using the
    Louvain algorithm. Proteins within the same community are more similar to
    each other than to proteins in other communities.

    Parameters
    ----------
    hits_tsv: Path
        Path to the input TSV file containing protein alignment hits.
        Expected columns: 'query', 'target', 'bits'.
    output_file : Path
        Path to the output TSV file where community membership will be saved.
        Will contain columns: 'protein_id', 'chg_id' (community/cluster ID).
    min_bitscore : int, optional
        Minimum bitscore threshold for including an edge in the network.
        Higher values create more stringent networks with stronger connections.
        Default is 100.
        - If clusters are too large, increase this value (e.g., 150)
        - If clusters are too small/fragmented, decrease this value (e.g., 80)

    Returns
    -------
    int
        Number of communities detected.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If no hits pass the bitscore filter.
    """
    # 1. Load the hits file
    logger.info(f"Loading alignment data from {hits_tsv}")
    try:
        hits_df = pd.read_csv(hits_tsv, sep="\t")
    except Exception as e:
        logger.error(f"Error reading file: {e}")
        raise FileNotFoundError(f"Could not read input file: {hits_tsv}") from e

    # 2. Filter hits to build the network
    logger.info(f"Original number of hits: {len(hits_df)}")
    filtered_df = hits_df[hits_df["bits"] > min_bitscore]
    logger.info(f"Hits after filtering (bitscore > {min_bitscore}): {len(filtered_df)}")

    if len(filtered_df) == 0:
        raise ValueError(
            f"No hits passed the filter (bitscore > {min_bitscore}). "
            f"Please lower the min_bitscore threshold."
        )

    # 3. Create the graph
    logger.info("Building the network graph...")
    # We tell networkx that 'query' is node 1, 'target' is node 2,
    # and the 'bits' column should be the 'weight' of the edge.
    G = nx.from_pandas_edgelist(filtered_df, source="query", target="target", edge_attr="bits")

    # 4. Run Louvain Community Detection
    logger.info("Running Louvain community detection algorithm...")
    # The Louvain algorithm uses the 'bits' weight to make decisions:
    # stronger connections (higher bitscore) are more likely
    # to be kept in the same community.
    communities = nx.community.louvain_communities(G, weight="bits")

    # 5. Save the results
    logger.info(f"Saving community membership to {chg_membership_tsv}")

    # Ensure output directory exists
    chg_membership_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Convert the list of sets into the two-column format
    # (protein_id, chg_id) where chg_id is the Community Homology Group ID
    with open(chg_membership_tsv, "w") as f:
        f.write("protein_id\tchg_id\n")

        # 'enumerate' gives us a unique chg_id (0, 1, 2, ...) for each community
        for chg_id, community_set in enumerate(communities):
            # 'community_set' is a set like {'proteinA', 'proteinB', ...}
            for protein in community_set:
                f.write(f"{protein}\t{chg_id}\n")

    num_communities = len(communities)
    logger.info(f"Community detection complete. Found {num_communities} communities.")
