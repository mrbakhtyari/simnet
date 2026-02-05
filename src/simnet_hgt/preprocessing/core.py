"""
I/O utilities for FASTA files.

This module handles merging protein sequences from multiple genome files
with systematic naming and metadata tracking.
"""

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

logger = logging.getLogger(__name__)


@dataclass
class ProteinMetadata:
    """Metadata for a protein sequence."""

    protein_id: str
    genome_accession: str
    taxon: str


def sanitize_text(text: str) -> str:
    """
    Sanitize text by replacing problematic characters.

    Args:
        text: Original text

    Returns:
        Sanitized text with spaces, pipes, commas, and other special chars replaced
    """
    return (
        text.replace(" ", "_")
        .replace("|", "_")
        .replace(":", "_")
        .replace(";", "_")
        .replace(",", "_")
        .replace("(", "_")
        .replace(")", "_")
    )


def extract_genome_accession(filename: Path) -> str:
    """
    Extract genome accession from protein filename.

    Args:
        filename: Path to protein file (e.g., NC_000854.2_proteins.faa)

    Returns:
        Genome accession (e.g., NC_000854.2)
    """
    return filename.stem.replace("_proteins", "")


def build_protein_id(
    taxon: str,
    genome_accession: str,
    protein_id: str,
    product: str,
    index: int,
) -> str:
    """
    Build standardized protein identifier.

    Format: TAXON_GenomeID_ProteinID_Product_Index

    Args:
        taxon: Taxonomic category (e.g., archaea, bacteria)
        genome_accession: Genome accession number
        protein_id: Original protein ID
        product: Protein product description
        index: Index counter for this genome file

    Returns:
        Formatted protein identifier
    """
    # Sanitize all components
    taxon_clean = sanitize_text(taxon)
    genome_clean = sanitize_text(genome_accession)
    protein_clean = sanitize_text(protein_id)
    product_clean = sanitize_text(product)

    return f"{taxon_clean}_{genome_clean}_{protein_clean}_{product_clean}_{index}"


def extract_product_from_description(description: str) -> str:
    """
    Extract product name from FASTA description.

    Args:
        description: Full FASTA description line

    Returns:
        Product description or 'hypothetical_protein' if not found
    """
    # Description format: "ID product_description"
    # Split by first space and take everything after
    parts = description.split(maxsplit=1)
    if len(parts) > 1:
        return parts[1]
    return "hypothetical_protein"


def collect_protein_files(
    raw_genomes_dir: Path,
    categories: list[str],
) -> list[tuple[str, Path]]:
    """
    Collect all protein files from specified categories.

    Args:
        raw_genomes_dir: Root directory containing genome data
        categories: List of category subdirectories to process

    Returns:
        List of (category, protein_file_path) tuples
    """
    protein_files = []

    for category in categories:
        category_dir = raw_genomes_dir / category
        if not category_dir.exists():
            logger.warning(f"Category directory not found: {category_dir}")
            continue

        for faa_file in category_dir.glob("*_proteins.faa"):
            protein_files.append((category, faa_file))

    return protein_files


def process_genome_proteins(
    taxon: str,
    protein_file: Path,
    fasta_writer,
    csv_writer,
) -> int:
    """
    Process all proteins from a single genome file.

    Args:
        taxon: Taxonomic category
        protein_file: Path to protein FASTA file
        fasta_writer: File handle for writing FASTA sequences
        csv_writer: CSV writer for writing metadata

    Returns:
        Number of proteins processed
    """
    genome_accession = extract_genome_accession(protein_file)
    protein_count = 0

    try:
        for record in SeqIO.parse(protein_file, "fasta"):
            # Increment counter for this genome (1-indexed)
            protein_count += 1

            # Extract components
            original_id = record.id
            product = extract_product_from_description(record.description)

            # Build new protein ID
            new_protein_id = build_protein_id(
                taxon=taxon,
                genome_accession=genome_accession,
                protein_id=original_id,
                product=product,
                index=protein_count,
            )

            # Write to FASTA
            fasta_writer.write(f">{new_protein_id}\n{record.seq!s}\n")

            # Write metadata to CSV
            csv_writer.writerow([new_protein_id, genome_accession, taxon])

    except Exception as e:
        logger.error(f"Error processing {protein_file}: {e}")

    return protein_count


def merge_and_map_proteins(
    raw_genomes_dir: Path,
    output_dir: Path,
    categories: list[str],
) -> tuple[int, int]:
    """
    Merge protein sequences and create metadata mapping in a single pass.

    This function processes all genome protein files, creating:
    1. all_proteins.faa - Merged FASTA with standardized IDs
    2. protein_metadata.csv - Mapping of protein_id -> genome_accession, taxon

    Protein IDs follow format: TAXON_GenomeID_ProteinID_Product_Index

    Args:
        raw_genomes_dir: Root directory containing genome category folders
        output_dir: Directory for output files
        categories: List of category names to process

    Returns:
        Tuple of (total_proteins, total_genomes) processed
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_output = output_dir / "all_proteins.faa"
    csv_output = output_dir / "protein_metadata.csv"

    logger.info("Starting protein merge and mapping")
    logger.info(f"Input directory: {raw_genomes_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Categories: {', '.join(categories)}")

    # Collect all protein files
    protein_files = collect_protein_files(raw_genomes_dir, categories)
    total_genomes = len(protein_files)
    total_proteins = 0

    logger.info(f"Found {total_genomes:,} genome files to process")

    # Open output files
    with (
        open(fasta_output, "w", encoding="utf-8") as fasta_out,
        open(csv_output, "w", encoding="utf-8", newline="") as csv_out,
    ):
        # Initialize CSV writer
        csv_writer = csv.writer(csv_out)
        csv_writer.writerow(["protein_id", "genome_accession", "taxon"])

        # Process each genome file
        show_progress = logger.isEnabledFor(logging.INFO)
        for taxon, protein_file in tqdm(
            protein_files,
            desc="Processing genomes",
            disable=not show_progress,
        ):
            count = process_genome_proteins(
                taxon=taxon,
                protein_file=protein_file,
                fasta_writer=fasta_out,
                csv_writer=csv_writer,
            )
            total_proteins += count

    logger.info(f"Merged {total_proteins:,} proteins from {total_genomes:,} genomes")
    logger.info(f"FASTA output: {fasta_output}")
    logger.info(f"Metadata output: {csv_output}")

    return total_proteins, total_genomes
