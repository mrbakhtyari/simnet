import logging
from pathlib import Path

from simnet_hgt.network_analysis import reduce_bipartite_network


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root = Path(__file__).parents[1] / "datasets/test"
    print(f"Using data root directory: {root}")

    networks_dir = root / "04_networks"

    connected_components_tsv = networks_dir / "connected_components.tsv"

    nodes = reduce_bipartite_network(
        bipartite_nodes_csv=networks_dir / "bipartite_network_nodes.csv",
        bipartite_edges_csv=networks_dir / "bipartite_network_edges.csv",
        connected_components_tsv=connected_components_tsv,
    )

    nodes_file = networks_dir / "genome_nodes.txt"
    with open(nodes_file, "w") as f:
        for node in nodes:
            f.write(f"{node}\n")
    logging.info(f"Saved {len(nodes)} genome nodes to {nodes_file}")


if __name__ == "__main__":
    main()
