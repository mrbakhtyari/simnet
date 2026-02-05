import logging
from pathlib import Path

from simnet_hgt.network_analysis import build_bipartite_network


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = Path(__file__).parents[1] / "datasets/test"
    print(f"Using data root directory: {root}")

    preprocessing_dir = root / "02_preprocessing"
    networks_dir = root / "04_networks"

    # Step 5: Bipartite graph construction and analysis
    connected_components_tsv = networks_dir / "connected_components.tsv"
    protein_metadata_csv = preprocessing_dir / "protein_metadata.csv"

    build_bipartite_network(
        connected_components_tsv=connected_components_tsv,
        protein_metadata_csv=protein_metadata_csv,
        output_dir=networks_dir,
    )


if __name__ == "__main__":
    main()
