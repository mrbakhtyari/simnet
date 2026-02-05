"""
Preprocessing module for genome protein data.

This module handles the preprocessing of genome protein sequences,
including merging FASTA files and creating metadata mappings.
"""

from .core import merge_and_map_proteins
from .input_parsing import parse_inputfile, prepare_raw_genomes

__all__ = [
    "merge_and_map_proteins",
    "parse_inputfile",
    "prepare_raw_genomes",
]
