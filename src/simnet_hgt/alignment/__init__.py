"""MMseqs2 module for protein sequence alignment."""

from .search import SearchParams, align_sequences

__all__ = [
    "SearchParams",
    "align_sequences",
]
