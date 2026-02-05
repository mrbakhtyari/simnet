import logging
import os
from pathlib import Path

import dotenv

from .pairwise_distance import compute_pairwise_distances
from .tree_builder import build_species_tree

dotenv.load_dotenv()


def build_phylogeny_pipeline(
    genome_accessions: list[str],
    output_path: Path,
    prefix: str = "",
) -> None:
    """
    Build species tree and compute pairwise distances for a set of genomes.

    Args:
        genome_accessions: List of genome accession IDs
        output_dir: Directory where output files will be saved
        prefix: Prefix for output files (e.g., 'all' or 'reduced')
        email: NCBI email for API access
        api_key: NCBI API key

    Returns:
        Tuple of (mean internal nodes, leaf count)
    """
    output_path.mkdir(parents=True, exist_ok=True)

    # Define output paths

    # Build species tree
    logging.info(f"Building {prefix} species tree with {len(genome_accessions)} genomes...")
    newick_path = build_species_tree(
        genome_accessions=genome_accessions,
        output_path=output_path,
        prefix=prefix,
        api_key=os.getenv("NCBI_API_KEY"),
        email=os.getenv("NCBI_EMAIL"),
    )

    # # Compute pairwise distances
    logging.info(f"Computing pairwise distances for {prefix} tree...")
    compute_pairwise_distances(
        newick_path=newick_path,
        output_path=output_path,
        prefix=prefix,
    )
