"""Utilities for reading and writing network data (nodes/edges) from/to CSV."""

import csv
from collections.abc import Iterable
from pathlib import Path

from .network import NetworkEdge, NetworkNode

# Color scheme for different categories
CATEGORY_COLORS = {
    "bacteria": "#2ca02c",  # green
    "archaea": "#ffcc00",  # yellow
    "plasmid": "#9467bd",  # purple
    "virus": "#d62728",  # red
}


def read_nodes(path: Path) -> list[NetworkNode]:
    """Read node metadata from CSV."""
    nodes: list[NetworkNode] = []
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            node_id = row.get("node_id") or row.get("genome_accession")
            if not node_id:
                continue
            nodes.append(
                NetworkNode(
                    node_id=node_id,
                    node_type=row.get("type", "genome" if row.get("category") else "CHG"),
                    category=row.get("category", ""),
                    color=row.get("color", "#bbbbbb"),
                    twin_group=row.get("twin_group", ""),
                    is_articulation=row.get("is_articulation", "False") == "True",
                )
            )
    return nodes


def read_edges(path: Path) -> list[NetworkEdge]:
    """Read edges from CSV."""
    edges: list[NetworkEdge] = []
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            edges.append(
                NetworkEdge(
                    source=row["source"],
                    target=row["target"],
                    weight=int(row.get("weight", 1)),
                )
            )
    return edges


def iter_edges_simple(path: Path) -> Iterable[tuple[str, str, int]]:
    """Iterate edges as simple tuples (for memory efficiency)."""
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            yield row["source"], row["target"], int(row.get("weight", 1))


def write_nodes(
    path: Path, nodes: list[NetworkNode], include_reduction_fields: bool = False
) -> None:
    """Write nodes to CSV.

    Args:
        path: Output CSV path
        nodes: List of NetworkNode objects
        include_reduction_fields: If True, include twin_group and is_articulation columns
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        if include_reduction_fields:
            w.writerow(["node_id", "type", "category", "color", "twin_group", "is_articulation"])
            for node in nodes:
                w.writerow(
                    [
                        node.node_id,
                        node.node_type,
                        node.category,
                        node.color,
                        node.twin_group,
                        str(node.is_articulation),
                    ]
                )
        else:
            w.writerow(["node_id", "type", "category", "color"])
            for node in nodes:
                w.writerow([node.node_id, node.node_type, node.category, node.color])


def write_edges(path: Path, edges: list[NetworkEdge]) -> None:
    """Write edges to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "target", "weight"])
        for edge in edges:
            w.writerow([edge.source, edge.target, edge.weight])


def load_protein_metadata(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Read protein_metadata.csv and extract both mappings.

    Args:
        path: Path to protein_metadata.csv file

    Returns:
        Tuple of (protein_to_genome, genome_to_category) dictionaries
    """
    protein_to_genome: dict[str, str] = {}
    genome_to_taxon: dict[str, str] = {}

    with open(path, encoding="utf-8") as f:
        r = csv.reader(f)
        header = next(r, None)

        # Verify header or assume it's present
        if header and header[0] == "protein_id":
            # Expected format: protein_id, genome_accession, taxon
            for row in r:
                protein_id, genome_accession, taxon = row[0], row[1], row[2]
                protein_to_genome[protein_id] = genome_accession
                genome_to_taxon[genome_accession] = taxon.lower()
        else:
            # No header, rewind and parse
            f.seek(0)
            r = csv.reader(f)
            for row in r:
                protein_id, genome_accession, taxon = row[0], row[1], row[2]
                protein_to_genome[protein_id] = genome_accession
                genome_to_taxon[genome_accession] = taxon.lower()

    return protein_to_genome, genome_to_taxon


def iter_cluster_members(chg_membership_tsv: Path) -> Iterable[tuple[str, str]]:
    """Iterate (protein_id, chg_id) tuples from CHG membership TSV."""
    with open(chg_membership_tsv, encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r, None)
        if header and header[1] != "chg_id":
            # no header, rewind
            f.seek(0)
            r = csv.reader(f, delimiter="\t")
        yield from ((protein_id, chg_id) for protein_id, chg_id in r)


def iter_protein_edges(
    hits_tsv: Path, has_header: bool = True, max_lines: int | None = None
) -> Iterable[tuple[str, str]]:
    """Iterate protein-protein edges from TSV file.

    Args:
        hits_tsv: Path to TSV with [query, target, ...] columns
        has_header: Whether file has a header row
        max_lines: Optional limit on number of lines to read

    Yields:
        (query, target) tuples
    """
    MIN_COLS = 2
    with open(hits_tsv, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        first = True
        for i, row in enumerate(reader):
            if first and has_header:
                first = False
                continue
            if not row or len(row) < MIN_COLS:
                continue
            if max_lines is not None and i > max_lines:
                break
            yield row[0], row[1]


def get_category_color(category: str) -> str:
    """Get visualization color for a genomic category."""
    return CATEGORY_COLORS.get(category.lower(), "#7f7f7f")
