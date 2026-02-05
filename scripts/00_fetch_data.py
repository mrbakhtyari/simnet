import logging
from pathlib import Path

import dotenv

from simnet_hgt.data_fetching import fetch_genome_data

dotenv.load_dotenv()


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = "datasets"
    # Step 0: Fetch test data from NCBI
    references_dir = Path(root) / "00_references"
    raw_genomes_dir = Path(root) / "01_raw_genomes"

    fetch_genome_data(references_dir=references_dir, output_dir=raw_genomes_dir)


if __name__ == "__main__":
    main()
