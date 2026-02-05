"""
File storage for genome data.
"""

import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


class GenomeStorage:
    """Handles saving files to disk."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)

    def save_fasta(self, data: str, filename: str, category: str) -> Path:
        """Save FASTA file."""
        filepath = self._get_filepath(filename, category, ".fasta")
        filepath.write_text(data, encoding="utf-8")
        return filepath

    def save_genbank(self, data: str, filename: str, category: str) -> Path:
        """Save GenBank file."""
        filepath = self._get_filepath(filename, category, ".gb")
        filepath.write_text(data, encoding="utf-8")
        return filepath

    def save_proteins(self, proteins: list[str], filename: str, category: str) -> Path:
        """Save protein FASTA file."""
        filepath = self._get_filepath(filename, category, "_proteins.faa")
        filepath.write_text("".join(proteins), encoding="utf-8")
        return filepath

    def save_genes(self, genes: list[dict[str, Any]], filename: str, category: str) -> Path:
        """Save gene JSON file."""
        filepath = self._get_filepath(filename, category, "_genes.json")
        filepath.write_text(json.dumps(genes, indent=2), encoding="utf-8")
        return filepath

    def save_metadata(self, results: list[Any], category: str) -> Path:
        """Save metadata Excel file."""
        # Convert to DataFrame
        records = [r.to_dict() for r in results]
        df = pd.DataFrame(records)

        # Define expected column order based on category
        if category in ("archaea", "bacteria"):
            expected_cols = [
                "UniqID",
                "Organism/Name",
                "TaxID",
                "BioProject Accession",
                "BioProject ID",
                "Group",
                "SubGroup",
                "Size (Mb)",
                "GC%",
                "Chromosomes/RefSeq",
                "Chromosomes/INSDC",
                "Plasmids/RefSeq",
                "Plasmids/INSDC",
                "WGS",
                "Scaffolds",
                "Genes",
                "Proteins",
                "Release Date",
                "Modify Date",
                "Status",
                "Center",
                "BioSample Accession",
                "Assembly Accession",
                "Body Products",
                "Biotic Relationship",
                "Viral Subgroup",
                "Temperature Optimum",
                "Symbiotic Physical Interaction",
                "Symbiont Name",
                "Cell Length",
                "Motility",
                "Commercial Strain Comments",
                "Gram Stain",
                "Ph",
                "Oxygen Requirement",
                "Metabolism",
                "Diseases",
                "Commercial Strain",
                "Phenotypes",
                "Cell Diameter",
                "Cell Arrangements",
                "Temperature Range",
                "Sporulation",
                "Energy Sources",
                "Salinity",
                "Habitats",
                "Symbiont Taxon ID",
                "Color",
                "Symbiotic Relationship",
                "Viral Group",
                "Cell Shape",
                "Kingdom",
                "Phylum",
            ]
        elif category == "plasmid":
            expected_cols = [
                "UniqID",
                "Organism/Name",
                "TaxID",
                "RefSeq",
                "GI",
                "GenBank",
                "HostName",
            ]
        elif category == "virus":
            expected_cols = [
                "UniqID",
                "Organism/Name",
                "RefSeq_Accession",  # Add RefSeq accession column
                "TaxID",
                "BioProject Accession",
                "BioProject ID",
                "Group",
                "SubGroup",
                "GC%",
                "Genes",
                "Proteins",
                "Release Date",
                "Modify Date",
                "Status",
                "Size (Kb)",
                "Host",
            ]
        else:
            expected_cols = []

        # Reorder columns: expected columns first, then any extras (prefixed with _)
        existing_expected = [col for col in expected_cols if col in df.columns]
        operational_cols = [col for col in df.columns if col.startswith("_")]
        other_cols = [
            col
            for col in df.columns
            if col not in existing_expected and col not in operational_cols
        ]

        final_columns = existing_expected + other_cols + operational_cols
        df = df[final_columns]

        # Save to category root (not subdirectory)
        filepath = self.output_dir / f"{category}_metadata.xlsx"
        df.to_excel(filepath, index=False)

        logger.info(f"Metadata saved: {filepath.name}")
        return filepath

    def _get_filepath(self, filename: str, category: str, suffix: str) -> Path:
        """Build filepath with safe filename."""
        # Don't sanitize - NCBI accessions are already filesystem-safe
        # Just ensure no directory traversal and limit length
        safe_name = filename.replace("/", "_").replace("\\", "_")[:200]
        category_dir = self.output_dir / category
        category_dir.mkdir(parents=True, exist_ok=True)
        return category_dir / f"{safe_name}{suffix}"

    @staticmethod
    def _sanitize_filename(name: str) -> str:
        """
        Convert to safe filename.

        NCBI accessions like NC_000852.5 or NT_167350.1 are already filesystem-safe.
        Only replace directory separators to prevent path traversal.
        """
        safe = name.replace("/", "_").replace("\\", "_")
        return safe[:200]
