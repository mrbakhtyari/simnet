"""
Dataset configuration for NCBI fetching.

Moved from ncbi_fetcher/config.py - NO LOGIC CHANGES.
"""

import logging
import re
from collections.abc import Hashable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from .genome import OrganismRecord

logger = logging.getLogger(__name__)


@dataclass
class DatasetConfig:
    """Configuration for a dataset to fetch from NCBI."""

    # Core fields
    name: str
    excel_file: Path

    # Column mapping
    accession_column: str
    uniq_id_column: str = "No."
    organism_column: str = "Organism/Name"
    taxid_column: str = "TaxID"
    plasmid_column: str | None = None
    genes_column: str | None = None
    proteins_column: str | None = None

    # Behavior flags
    combine_chromosomes: bool = True
    include_plasmids: bool = False
    fetch_proteins: bool = True
    fetch_genes: bool = True
    validate_counts: bool = False

    # Direct accessions (alternative to Excel)
    _direct_accessions: list[str] | None = field(default=None, repr=False)

    # ...existing classmethods (for_archaea, for_bacteria, etc.)...

    @classmethod
    def for_archaea(cls, excel_file: Path) -> "DatasetConfig":
        """Create config for Archaea dataset (CorelS1.xlsx)."""
        return cls(
            name="archaea",
            excel_file=excel_file,
            accession_column="Chromosomes/RefSeq",
            uniq_id_column="UniqID",
            plasmid_column="Plasmids/RefSeq",
            include_plasmids=True,
            combine_chromosomes=True,
            genes_column="Genes",
            proteins_column="Proteins",
        )

    @classmethod
    def for_bacteria(cls, excel_file: Path) -> "DatasetConfig":
        """Create config for Bacteria dataset (CorelS2.xlsx)."""
        return cls(
            name="bacteria",
            excel_file=excel_file,
            accession_column="Chromosomes/RefSeq",
            uniq_id_column="UniqID",
            plasmid_column="Plasmids/RefSeq",
            include_plasmids=True,
            combine_chromosomes=True,
            genes_column="Genes",
            proteins_column="Proteins",
        )

    @classmethod
    def for_plasmid(cls, excel_file: Path) -> "DatasetConfig":
        """Create config for Plasmid dataset (CorelS3.xlsx)."""
        try:
            header_cols = set(pd.read_excel(excel_file, nrows=0).columns)
        except Exception:
            header_cols = set()

        return cls(
            name="plasmid",
            excel_file=excel_file,
            accession_column="RefSeq",
            uniq_id_column="UniqID",
            include_plasmids=False,
            combine_chromosomes=False,
            genes_column="Genes" if "Genes" in header_cols else None,
            proteins_column="Proteins" if "Proteins" in header_cols else None,
            validate_counts=bool("Genes" in header_cols or "Proteins" in header_cols),
        )

    @classmethod
    def for_virus(cls, excel_file: Path) -> "DatasetConfig":
        """Create config for Virus dataset (CorelS4.xlsx)."""
        return cls(
            name="virus",
            excel_file=excel_file,
            accession_column="BioProject Accession",
            uniq_id_column="UniqID",
            include_plasmids=False,
            combine_chromosomes=False,
            validate_counts=False,
        )

    @classmethod
    def from_accessions(  # noqa: PLR0913
        cls,
        name: str,
        accessions: list[str],
        combine_chromosomes: bool = False,
        include_plasmids: bool = False,
        fetch_proteins: bool = True,
        fetch_genes: bool = True,
    ) -> "DatasetConfig":
        """Create config from direct accession list."""
        config = cls(
            name=name,
            excel_file=Path(""),
            accession_column="",
            combine_chromosomes=combine_chromosomes,
            include_plasmids=include_plasmids,
            fetch_proteins=fetch_proteins,
            fetch_genes=fetch_genes,
            validate_counts=False,
        )
        config._direct_accessions = accessions
        return config

    def load_records(self) -> list[OrganismRecord]:
        """Load organism records from Excel file or direct accessions."""
        # Check if using direct accessions
        if self._direct_accessions is not None:
            return self._load_records_from_accessions()

        # Otherwise use Excel file
        if not self.excel_file.exists():
            logger.error(f"Excel file not found: {self.excel_file}")
            return []

        logger.info(f"Loading records from {self.excel_file.name}")
        df = pd.read_excel(self.excel_file)

        records = []
        for idx, row in df.iterrows():
            try:
                record = self._create_record_from_row(row, idx)
                if record and record.accessions:
                    records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse row {idx}: {e}")
                continue

        logger.info(f"Loaded {len(records)} records from {self.name} dataset")
        return records

    def _load_records_from_accessions(self) -> list[OrganismRecord]:
        """Load records from direct accession list."""
        if not self._direct_accessions:
            return []

        records = []
        for i, accession in enumerate(self._direct_accessions, start=1):
            record = OrganismRecord(
                uniq_id=f"{self.name}_{i}",
                organism="Unknown",
                taxid=None,
                category=self.name,
                accessions=[accession],
                additional_fields={},
            )
            records.append(record)

        return records

    @staticmethod
    def _split_accession_list(text: str) -> list[str]:
        """Split accession list by common delimiters."""
        if not text or pd.isna(text):
            return []

        # Split on commas, semicolons, pipes, or whitespace
        parts = re.split(r"[,;\s|]+", str(text).strip())

        # Filter out empty strings and placeholder values
        accessions = []
        for part in parts:
            cleaned = part.strip()
            if cleaned and cleaned not in {"-", "nan", "N/A", ""}:
                accessions.append(cleaned)

        return accessions

    def _extract_accessions(self, row: pd.Series) -> list[str]:
        """Extract accession(s) from row, optionally including plasmids."""
        accessions = []

        # Extract main accessions (chromosomes)
        acc_value = row.get(self.accession_column)
        if pd.notna(acc_value):
            accessions.extend(self._split_accession_list(str(acc_value)))

        # Optionally include plasmids
        if self.include_plasmids and self.plasmid_column:
            plasmid_value = row.get(self.plasmid_column)
            if pd.notna(plasmid_value):
                accessions.extend(self._split_accession_list(str(plasmid_value)))

        return accessions

    @staticmethod
    def _safe_int(value: Any) -> int | None:
        """
        Safely convert value to integer.

        Handles:
        - None/NaN -> None
        - Floats: 1234.0 -> 1234
        - Strings with commas: "1,234" -> 1234
        - Placeholder values like "-" -> None
        - Invalid values -> None
        """
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None

        s = str(value).strip()
        if not s or s in {"-", "nan", "N/A", ""}:
            return None

        # Remove all non-digit characters
        digits = re.sub(r"[^\d]", "", s)

        try:
            return int(digits) if digits else None
        except ValueError:
            return None

    @staticmethod
    def _normalize_taxid(value: Any) -> str | None:
        """
        Normalize TaxID to clean string.

        Excel may store as: 511145, 511145.0, "511145", "-", etc.
        Returns clean digit string or None.
        """
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None

        s = str(value).strip()
        if not s or s in {"-", "nan", "N/A", ""}:
            return None

        # Keep only digits
        digits = re.sub(r"[^\d]", "", s)
        return digits if digits else None

    def _create_record_from_row(self, row: pd.Series, idx: Hashable) -> OrganismRecord | None:
        """Create an OrganismRecord from a DataFrame row."""
        # Extract accessions first
        accessions = self._extract_accessions(row)
        if not accessions:
            logger.debug(f"Row {idx}: No valid accessions found")
            return None

        # Extract fields with safe conversions
        uniq_id = str(row.get(self.uniq_id_column, f"{self.name}_{idx}"))
        organism = str(row.get(self.organism_column, "Unknown"))
        taxid = self._normalize_taxid(row.get(self.taxid_column))

        # Extract metadata
        additional_fields = self._extract_additional_fields(row)

        # Preserve original accession columns for archaea/bacteria
        if self.name in ("archaea", "bacteria"):
            additional_fields["_original_chromosomes_refseq"] = row.get(self.accession_column, "")
            additional_fields["_original_chromosomes_insdc"] = row.get("Chromosomes/INSDC", "")
            if self.plasmid_column:
                additional_fields["_original_plasmids_refseq"] = row.get(self.plasmid_column, "-")
                additional_fields["_original_plasmids_insdc"] = row.get("Plasmids/INSDC", "-")

        return OrganismRecord(
            uniq_id=uniq_id,
            organism=organism,
            taxid=taxid,
            category=self.name,
            accessions=accessions,
            additional_fields=additional_fields,
        )

    def _extract_additional_fields(self, row: pd.Series) -> dict[str, Any]:
        """
        Extract minimal fields from Excel - only expected counts for validation.

        All other metadata will be fetched from NCBI.
        These are stored with _ prefix to distinguish from actual metadata.
        """
        additional_fields = {}
        row_dict = row.to_dict()

        # Store expected counts for validation (safely) with _ prefix
        if self.genes_column and self.genes_column in row_dict:
            genes_val = self._safe_int(row_dict[self.genes_column])
            if genes_val is not None:
                additional_fields["_expected_genes"] = genes_val

        if self.proteins_column and self.proteins_column in row_dict:
            proteins_val = self._safe_int(row_dict[self.proteins_column])
            if proteins_val is not None:
                additional_fields["_expected_proteins"] = proteins_val

        if self.name == "virus" and "BioProject ID" in row_dict:
            bpid_val = self._safe_int(row_dict["BioProject ID"])
            if bpid_val is not None:
                additional_fields["_excel_bioproject_id"] = str(bpid_val)

        return additional_fields
