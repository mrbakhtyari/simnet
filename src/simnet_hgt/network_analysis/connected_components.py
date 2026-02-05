"""Domain models and algorithms for protein clustering (CHG)."""

import csv
import logging
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .io import iter_protein_edges

logger = logging.getLogger(__name__)


@dataclass
class CHG:
    """Cluster of Homologous Genes - a group of similar proteins."""

    chg_id: str
    protein_ids: set[str] = field(default_factory=set)

    def add_protein(self, protein_id: str) -> None:
        """Add a protein to this cluster."""
        self.protein_ids.add(protein_id)

    @property
    def size(self) -> int:
        """Number of proteins in this cluster."""
        return len(self.protein_ids)


class UnionFind:
    """Union-Find data structure for computing connected components.

    This is a domain algorithm for clustering proteins based on homology edges.
    It's pure domain logic with no external dependencies.
    """

    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def find(self, x: str) -> str:
        """Find the root of the set containing x (with path compression)."""
        if self.parent.get(x, x) != x:
            self.parent[x] = self.find(self.parent[x])
        else:
            # Initialize unseen node
            self.parent.setdefault(x, x)
            self.rank.setdefault(x, 0)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        """Merge the sets containing a and b (with union by rank)."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        # Union by rank
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

    def get_all_nodes(self) -> set[str]:
        """Return all nodes that have been added to the structure."""
        return set(self.parent.keys())

    def get_components(self) -> dict[str, list[str]]:
        """Return all connected components as a dict of root -> members."""
        by_root: dict[str, list[str]] = defaultdict(list)
        for node in self.parent:
            root = self.find(node)
            by_root[root].append(node)
        return dict(by_root)


def detect_connected_components(
    hits_tsv: Path,
    connected_components_tsv: Path,
    has_header: bool = True,
) -> float:
    """Compute CHG membership via connected components.

    Args:
        hits_tsv: TSV with columns [query, target, ...]
        output_membership_tsv: Output TSV path with columns [protein_id, chg_id]
        has_header: Whether hits_tsv has a header row
    """
    connected_components_tsv.parent.mkdir(parents=True, exist_ok=True)

    # 1) Build Union-Find structure from protein edges
    # Load all edges first to isolate I/O vs computation
    edges = list(iter_protein_edges(hits_tsv, has_header))

    start_time = time.perf_counter()
    uf = UnionFind()
    for protein_a, protein_b in edges:
        uf.union(protein_a, protein_b)

    # 2) Get connected components
    components = uf.get_components()

    # 3) Create CHG objects and assign IDs (filter singleton clusters)
    chgs: list[CHG] = []
    for chg_idx, (_root, members) in enumerate(components.items(), start=1):
        chg = CHG(chg_id=f"CHG_{chg_idx}")
        for protein_id in members:
            chg.add_protein(protein_id)
        chgs.append(chg)
    end_time = time.perf_counter()

    # 4) Write membership file
    num_proteins = 0
    with open(connected_components_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["protein_id", "chg_id"])
        for chg in chgs:
            for protein_id in chg.protein_ids:
                writer.writerow([protein_id, chg.chg_id])
                num_proteins += 1
    logger.info(f"Computed {len(chgs)} non-trivial CHGs with total {num_proteins} proteins. ")
    return end_time - start_time
