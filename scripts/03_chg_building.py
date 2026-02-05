import logging
from pathlib import Path

from simnet_hgt.network_analysis import detect_connected_components


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = Path(__file__).parents[1] / "datasets/test"
    print(f"Using data root directory: {root}")

    alignment_dir = root / "03_alignment"
    networks_dir = root / "04_networks"

    # Step 4: CHG building
    hits_tsv = alignment_dir / "hits_cleaned.tsv"
    connected_components_tsv = networks_dir / "connected_components.tsv"
    detect_connected_components(
        hits_tsv=hits_tsv,
        connected_components_tsv=connected_components_tsv,
    )


if __name__ == "__main__":
    main()
