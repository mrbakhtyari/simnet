"""HGT scoring script.

Computes HGT (Horizontal Gene Transfer) scores from protein similarity data
and performs statistical analysis to identify high-confidence candidates.
"""

import json
import logging
from pathlib import Path

from simnet_hgt.hgt import (
    HgtParams,
    ScoringMethod,
    analyze_hgt_scores,
    compute_hgt_scores,
    compute_similarity_mean,
)

SCORING_METHOD = ScoringMethod.EXPONENTIAL


def main() -> None:
    """Run the HGT scoring and analysis pipeline."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)

    root = Path(__file__).parents[1] / "datasets/test"
    logger.info("Using data root directory: %s", root)

    preprocessing_dir = root / "02_preprocessing"
    alignment_dir = root / "03_alignment"
    networks_dir = root / "04_networks"
    phylogeny_dir = root / "05_phylogeny"
    hgt_output_dir = root / "06_hgt"

    hits_tsv = alignment_dir / "hits_cleaned.tsv"

    reduced_chg_membership_tsv = networks_dir / "reduced_chg_membership.tsv"
    pairwise_distances_tsv = phylogeny_dir / "reduced_pairwise_distances.tsv"
    reduced_phylogeny_stats = phylogeny_dir / "reduced_phylogeny_stats.json"
    protein_metadata_csv = preprocessing_dir / "protein_metadata.csv"
    genome_to_taxid = phylogeny_dir / "genome_to_organism_name.csv"

    with reduced_phylogeny_stats.open() as f:
        data = json.load(f)
    n_int_mean = data.get("n_int_mean", 0.0)

    similarity_mean_json = alignment_dir / "similarity_mean.json"
    hits_mean = compute_similarity_mean(hits_file=hits_tsv, cache_file=similarity_mean_json)

    params = HgtParams(
        n_int_mean=n_int_mean,
        hits_mean=hits_mean,
        min_phylogenetic_distance=2,
        scoring_method=SCORING_METHOD,
    )

    compute_hgt_scores(
        hits_tsv=hits_tsv,
        file_paths={
            "protein_to_genome": protein_metadata_csv,
            "protein_to_chg": reduced_chg_membership_tsv,
            "phylogeny_distances": pairwise_distances_tsv,
            "genome_to_taxid": genome_to_taxid,
            "output": hgt_output_dir,
        },
        params=params,
    )

    analyze_hgt_scores(hgt_output_dir / "hgt_scores.tsv", hgt_output_dir)


if __name__ == "__main__":
    main()
