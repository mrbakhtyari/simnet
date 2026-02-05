"""Domain models for network entities."""

from dataclasses import dataclass, field


@dataclass
class NetworkNode:
    """Represents a node in a biological network."""

    node_id: str
    node_type: str  # "genome" or "CHG"
    category: str = ""  # archaea, bacteria, plasmid, virus (for genomes only)
    color: str = "#bbbbbb"  # Visualization color
    twin_group: str = ""  # For reduced graphs
    is_articulation: bool = False  # Public goods candidate


@dataclass
class NetworkEdge:
    """Represents an edge between two nodes."""

    source: str
    target: str
    weight: int = 1


@dataclass
class BipartiteGraph:
    """Represents a genome-CHG bipartite network."""

    nodes: list[NetworkNode] = field(default_factory=list)
    edges: list[NetworkEdge] = field(default_factory=list)

    @property
    def genome_nodes(self) -> list[NetworkNode]:
        """Return only genome nodes."""
        return [n for n in self.nodes if n.node_type == "genome"]

    @property
    def chg_nodes(self) -> list[NetworkNode]:
        """Return only CHG nodes."""
        return [n for n in self.nodes if n.node_type == "CHG"]

    @property
    def num_nodes(self) -> int:
        """Total number of nodes."""
        return len(self.nodes)

    @property
    def num_edges(self) -> int:
        """Total number of edges."""
        return len(self.edges)
