"""Data loading functions for HGT analysis.

This module provides functions to load various data files required for
horizontal gene transfer (HGT) analysis, including phylogenetic distances,
protein-to-genome mappings, and HGT scores.
"""

import csv
import logging
from pathlib import Path

import pandas as pd

from .models import (
    GenomeAccession,
    GenomeToTaxId,
    PhylogenyDistances,
    ProteinToChg,
    ProteinToGenome,
)

logger = logging.getLogger(__name__)


def load_phylogeny_distances(phylogeny_distances_file: Path) -> PhylogenyDistances:
    """Load phylogenetic distance (N_int) lookup table.

    Args:
        phylogeny_distances_file: Path to TSV file with columns:
            genome_a, genome_b, n_int

    Returns:
        Dictionary mapping (genomeA, genomeB) tuples to N_int values.
        Keys are canonically ordered (alphabetically) to ensure consistent lookup.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    logger.debug("Loading phylogeny distances from %s", phylogeny_distances_file)
    phylogeny_distances: PhylogenyDistances = {}

    with phylogeny_distances_file.open(encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)  # Skip header
        for row in reader:
            if len(row) >= 3:  # noqa: PLR2004
                genome_a, genome_b, n_int = row[0], row[1], int(row[2])
                # Canonical key ordering for consistent lookup
                key = (genome_a, genome_b) if genome_a < genome_b else (genome_b, genome_a)
                phylogeny_distances[key] = n_int

    logger.info("Loaded %d phylogeny distance pairs", len(phylogeny_distances))
    return phylogeny_distances


def load_protein_to_genome(protein_to_genome_file: Path) -> ProteinToGenome:
    """Load protein to genome mapping from file.

    Args:
        protein_to_genome_file: Path to CSV file with columns:
            protein_id, genome_accession

    Returns:
        Dictionary mapping protein IDs to genome accessions.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    logger.debug("Loading protein to genome mapping from %s", protein_to_genome_file)
    protein_to_genome: ProteinToGenome = {}

    with protein_to_genome_file.open(encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)  # Skip header
        for row in reader:
            if len(row) >= 2:  # noqa: PLR2004
                protein_to_genome[row[0]] = row[1]

    logger.info("Loaded %d protein-to-genome mappings", len(protein_to_genome))
    return protein_to_genome


def load_protein_to_chg(protein_to_chg_file: Path) -> ProteinToChg:
    """Load protein to CHG (Cluster of Homologous Genes) mapping from file.

    Args:
        protein_to_chg_file: Path to TSV file with columns:
            protein_id, chg_id

    Returns:
        Dictionary mapping protein IDs to CHG IDs.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    logger.debug("Loading protein to CHG mapping from %s", protein_to_chg_file)
    protein_to_chg: ProteinToChg = {}

    with protein_to_chg_file.open(encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)  # Skip header
        for row in reader:
            if len(row) >= 2:  # noqa: PLR2004
                protein_to_chg[row[0]] = row[1]

    logger.info("Loaded %d protein-to-CHG mappings", len(protein_to_chg))
    return protein_to_chg


def load_genome_to_taxid(genome_to_taxid_file: Path) -> GenomeToTaxId:
    """Load genome accession to taxid mapping from file.

    Args:
        genome_to_taxid_file: Path to CSV file with columns:
            genome_accession, organism_name, taxid

    Returns:
        Dictionary mapping genome accessions to taxids.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    logger.debug("Loading genome to taxid mapping from %s", genome_to_taxid_file)
    genome_to_taxid: GenomeToTaxId = {}

    with genome_to_taxid_file.open(encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)  # Skip header
        for row in reader:
            if len(row) >= 3:  # noqa: PLR2004
                # Columns: genome_accession, organism_name, taxid
                genome_to_taxid[row[0]] = row[2]

    logger.info("Loaded %d genome-to-taxid mappings", len(genome_to_taxid))
    return genome_to_taxid


def load_hgt_scores(hgt_scores_file: Path) -> pd.DataFrame:
    """Load computed HGT scores from a TSV file.

    Args:
        hgt_scores_file: Path to TSV file with HGT scores.

    Returns:
        DataFrame containing HGT scores with columns:
        protein_A, protein_B, CHG, TaxID_A, TaxID_B, N_int, s_AB, HGT_score

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    logger.debug("Loading HGT scores from %s", hgt_scores_file)
    df = pd.read_csv(hgt_scores_file, sep="\t")
    logger.info("Loaded %d HGT score records", len(df))
    return df


def get_canonical_genome_key(
    genome_a: GenomeAccession, genome_b: GenomeAccession
) -> tuple[GenomeAccession, GenomeAccession]:
    """Create a canonical key for genome pair lookup.

    Ensures consistent ordering for dictionary lookups regardless of
    which genome is specified first.

    Args:
        genome_a: First genome accession.
        genome_b: Second genome accession.

    Returns:
        Tuple of genome accessions in canonical (alphabetical) order.
    """
    return (genome_a, genome_b) if genome_a < genome_b else (genome_b, genome_a)
