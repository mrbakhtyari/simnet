import os
from pathlib import Path

from simnet_hgt.data_fetching import DatasetConfig, GenomeFetcher


def fetch_genome_data(references_dir: Path, output_dir: Path) -> None:
    """Fetch genome data from NCBI based on predefined datasets."""
    # # Initialize genome fetcher
    genome_fetcher = GenomeFetcher(
        output_dir=output_dir,
        email=os.getenv("NCBI_EMAIL"),
        api_key=os.getenv("NCBI_API_KEY"),
    )

    # Fetch records from each category of dataset
    archaea_config = DatasetConfig.for_archaea(references_dir / "CorelS1.xlsx")
    genome_fetcher.fetch_dataset(archaea_config)

    bacteria_config = DatasetConfig.for_bacteria(references_dir / "CorelS2.xlsx")
    genome_fetcher.fetch_dataset(bacteria_config)

    plasmid_config = DatasetConfig.for_plasmid(references_dir / "CorelS3.xlsx")
    genome_fetcher.fetch_dataset(plasmid_config)

    virus_config = DatasetConfig.for_virus(references_dir / "CorelS4.xlsx")
    genome_fetcher.fetch_dataset(virus_config)
