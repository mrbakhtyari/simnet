import logging
from pathlib import Path

from simnet_hgt.phylogeny import build_phylogeny_pipeline


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = Path(__file__).parents[1] / "datasets/test"
    print(f"Using data root directory: {root}")

    networks_dir = root / "04_networks"
    phylogeny_dir = root / "05_phylogeny"

    nodes_file = networks_dir / "genome_nodes.txt"
    if not nodes_file.exists():
        logging.error(
            f"Nodes file not found: {nodes_file}. Please run 06_network_reduction.py first."
        )
        return

    with open(nodes_file) as f:
        nodes = [line.strip() for line in f if line.strip()]

    build_phylogeny_pipeline(
        genome_accessions=nodes,
        output_path=phylogeny_dir,
        prefix="reduced",
    )


if __name__ == "__main__":
    main()
