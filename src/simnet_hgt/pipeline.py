import json
import logging
import time
from pathlib import Path

from simnet_hgt.alignment import SearchParams, align_sequences
from simnet_hgt.hgt import (
    HgtParams,
    ScoringMethod,
    analyze_hgt_scores,
    compute_hgt_scores,
    compute_similarity_mean,
    restore_original_names,
)
from simnet_hgt.network_analysis import (
    build_bipartite_network,
    detect_connected_components,
    reduce_bipartite_network,
)
from simnet_hgt.phylogeny import compute_pairwise_distances, prepare_phylogeny
from simnet_hgt.preprocessing import merge_and_map_proteins, prepare_raw_genomes

logger = logging.getLogger(__name__)


class HGTPipeline:
    """
    Orchestrates the HGT detection pipeline.
    """

    def __init__(
        self,
        dataset_path: Path,
        min_coverage: float = 0.8,
        min_identity: float = 0.59,
    ) -> None:
        self.dataset_path = dataset_path
        self.min_coverage = min_coverage
        self.min_identity = min_identity

        # Define directory structure
        self.dirs = {
            "input": dataset_path / "00_input",
            "raw_genomes": dataset_path / "01_raw_genomes",
            "preprocessing": dataset_path / "02_preprocessing",
            "alignment": dataset_path / "03_alignment",
            "networks": dataset_path / "04_networks",
            "phylogeny": dataset_path / "05_phylogeny",
            "hgt": dataset_path / "06_hgt",
        }
        self.timings: dict[str, float] = {}

    def run(self) -> None:
        """Execute the full pipeline."""
        logger.info(f"Starting HGT pipeline on {self.dataset_path}")
        start_time = time.perf_counter()

        name_mapping = self._step_0_prepare_inputs()
        self._step_1_preprocessing()
        self._step_2_alignment()

        # Core computation block
        self._step_3_chg_building()
        self._step_4_bipartite_analysis()
        self._step_5_pairwise_distances()
        self._step_6_hgt_scoring(name_mapping)

        self._step_final_postprocessing(name_mapping)
        total_elapsed = time.perf_counter() - start_time
        logger.info(f"Pipeline completed in {total_elapsed:.2f} seconds")
        self._print_timings()

    def _print_timings(self) -> None:
        """Print table of computation times."""
        print("\n" + "=" * 40)
        print(f"{'Step':<25} | {'Time (ms)':<10}")
        print("-" * 40)
        total_comp = 0.0
        for step, duration in self.timings.items():
            duration_ms = duration * 1000
            print(f"{step:<25} | {duration_ms:<10.4f}")
            total_comp += duration

        total_comp_ms = total_comp * 1000
        print("-" * 40)
        print(f"{'Total Computation':<25} | {total_comp_ms:<10.4f}")
        print("=" * 40 + "\n")

    def _step_0_prepare_inputs(self) -> dict[str, str]:
        """Prepare raw genomes and phylogeny input."""
        category = "unclassified"
        name_mapping = prepare_raw_genomes(
            input_dir=self.dirs["input"],
            raw_genomes_dir=self.dirs["raw_genomes"],
            category=category,
        )
        prepare_phylogeny(
            input_dir=self.dirs["input"],
            phylogeny_dir=self.dirs["phylogeny"],
            name_mapping=name_mapping,
        )
        return name_mapping

    def _step_1_preprocessing(self) -> None:
        """Merge and map proteins."""
        merge_and_map_proteins(
            raw_genomes_dir=self.dirs["raw_genomes"],
            output_dir=self.dirs["preprocessing"],
            categories=["unclassified"],
        )

    def _step_2_alignment(self) -> None:
        """Run sequence alignment."""
        params = SearchParams(
            min_coverage=self.min_coverage,
            min_identity=self.min_identity,
            coverage_mode=0,
            e_value=1e-5,
        )
        _, duration = align_sequences(
            fasta_path=self.dirs["preprocessing"] / "all_proteins.faa",
            output_dir=self.dirs["alignment"],
            params=params,
        )
        self.timings["Alignment (MMseqs2)"] = duration

    def _step_3_chg_building(self) -> None:
        """Detect connected components."""
        duration = detect_connected_components(
            hits_tsv=self.dirs["alignment"] / "hits_cleaned.tsv",
            connected_components_tsv=self.dirs["networks"] / "connected_components.tsv",
        )
        self.timings["CHG Detection"] = duration

    def _step_4_bipartite_analysis(self) -> None:
        """Build and reduce bipartite network."""
        connected_components_tsv = self.dirs["networks"] / "connected_components.tsv"
        protein_metadata_csv = self.dirs["preprocessing"] / "protein_metadata.csv"

        t1 = build_bipartite_network(
            connected_components_tsv=connected_components_tsv,
            protein_metadata_csv=protein_metadata_csv,
            output_dir=self.dirs["networks"],
        )

        nodes, t2 = reduce_bipartite_network(
            bipartite_nodes_csv=self.dirs["networks"] / "bipartite_network_nodes.csv",
            bipartite_edges_csv=self.dirs["networks"] / "bipartite_network_edges.csv",
            connected_components_tsv=connected_components_tsv,
        )
        self.timings["Network Construction"] = t1
        self.timings["Network Reduction"] = t2

        # Save nodes for debugging/reference
        nodes_file = self.dirs["networks"] / "genome_nodes.txt"
        with open(nodes_file, "w") as f:
            for node in nodes:
                f.write(f"{node}\n")

    def _step_5_pairwise_distances(self) -> None:
        """Compute pairwise phylogenetic distances."""
        duration = compute_pairwise_distances(
            newick_path=self.dirs["phylogeny"] / "tree.newick",
            output_path=self.dirs["phylogeny"],
            prefix="reduced",
        )
        self.timings["Phylogeny Distances"] = duration

    def _step_6_hgt_scoring(self, name_mapping: dict[str, str]) -> None:
        """Compute HGT scores."""
        # Load mean interaction stats
        reduced_phylogeny_stats = self.dirs["phylogeny"] / "reduced_phylogeny_stats.json"
        with reduced_phylogeny_stats.open() as f:
            data = json.load(f)
        n_int_mean = data.get("n_int_mean", 0.0)

        hits_tsv = self.dirs["alignment"] / "hits_cleaned.tsv"
        similarity_mean_json = self.dirs["alignment"] / "similarity_mean.json"
        hits_mean = compute_similarity_mean(hits_file=hits_tsv, cache_file=similarity_mean_json)

        hgt_params = HgtParams(
            n_int_mean=n_int_mean,
            hits_mean=hits_mean,
            min_phylogenetic_distance=2,
            scoring_method=ScoringMethod.RANKING_TH,
            rank_threshold=0.5,
        )

        _, duration = compute_hgt_scores(
            hits_tsv=hits_tsv,
            file_paths={
                "protein_to_genome": self.dirs["preprocessing"] / "protein_metadata.csv",
                "protein_to_chg": self.dirs["networks"] / "reduced_chg_membership.tsv",
                "phylogeny_distances": self.dirs["phylogeny"] / "reduced_pairwise_distances.tsv",
                "genome_to_taxid": self.dirs["phylogeny"] / "genome_to_organism_name.csv",
                "output": self.dirs["hgt"],
            },
            params=hgt_params,
        )
        self.timings["HGT Scoring"] = duration

        analyze_hgt_scores(self.dirs["hgt"] / "hgt_scores.tsv", self.dirs["hgt"])

    def _step_final_postprocessing(self, name_mapping: dict[str, str]) -> None:
        """Restore original organism names in output."""
        restore_original_names(
            hgt_output_dir=self.dirs["hgt"],
            name_mapping=name_mapping,
        )
