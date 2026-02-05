import logging
from pathlib import Path

from simnet_hgt.alignment import SearchParams, align_sequences


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = Path(__file__).parents[1] / "datasets/test"
    print(f"Using data root directory: {root}")

    preprocessing_dir = root / "02_preprocessing"
    alignment_dir = root / "03_alignment"

    # Step 2: MMseqs2 search for sequence alignment
    all_proteins_faa = preprocessing_dir / "all_proteins.faa"
    params = SearchParams(min_coverage=0.8, min_identity=0.95, coverage_mode=0, e_value=1e-5)
    align_sequences(
        fasta_path=all_proteins_faa,
        output_dir=alignment_dir,
        params=params,
    )


if __name__ == "__main__":
    main()
