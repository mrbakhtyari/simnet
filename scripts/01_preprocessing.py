import logging
from pathlib import Path

from simnet_hgt.preprocessing import merge_and_map_proteins


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = Path(__file__).parents[1] / "datasets/test"
    print(f"Using data root directory: {root}")

    preprocessing_dir = root / "02_preprocessing"
    raw_genomes_dir = root / "01_raw_genomes"

    # Step 1: Preprocessing
    categories = ["archaea", "bacteria", "plasmid", "virus"]
    merge_and_map_proteins(
        raw_genomes_dir=raw_genomes_dir,
        output_dir=preprocessing_dir,
        categories=categories,
    )


if __name__ == "__main__":
    main()
