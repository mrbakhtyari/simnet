"""
Data models for genome fetching.

Moved from ncbi_fetcher/models.py - NO LOGIC CHANGES.
"""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class OrganismRecord:
    """Represents a single organism record from reference Excel."""

    uniq_id: str
    organism: str
    taxid: str | None
    category: str
    accessions: list[str]
    additional_fields: dict[str, Any] = field(default_factory=dict)

    # Fields populated during fetching
    resolved_accessions: list[str] = field(default_factory=list)

    def get_primary_accession(self) -> str:
        """Get the primary accession for fetching."""
        return self.accessions[0] if self.accessions else ""

    def has_multiple_chromosomes(self) -> bool:
        """Check if organism has multiple chromosomes."""
        return len(self.accessions) > 1


@dataclass
class FetchResult:
    """Result of fetching sequence data for an organism."""

    record: OrganismRecord
    status: str  # 'success' or 'failed'
    error: str | None = None

    # Sequence information
    sequence_length: int = 0
    num_genes: int = 0
    num_proteins: int = 0

    # File paths
    fasta_file: str | None = None
    genbank_file: str | None = None
    protein_file: str | None = None
    gene_file: str | None = None

    def is_success(self) -> bool:
        """Check if fetch was successful."""
        return self.status in ("success", "cached")

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for metadata export."""
        # Start with all additional fields (which now include all GenBank metadata)
        result = dict(self.record.additional_fields)

        # Override/add key fields to ensure they're present
        result["UniqID"] = self.record.uniq_id

        # Only add these if not already in metadata
        if "Organism/Name" not in result:
            result["Organism/Name"] = self.record.organism
        if "TaxID" not in result:
            result["TaxID"] = self.record.taxid

        # Add operational fields (not from reference Excel)
        result["_category"] = self.record.category
        result["_original_accessions"] = ",".join(self.record.accessions)
        result["_status"] = self.status
        result["_error"] = self.error
        result["_sequence_length"] = self.sequence_length
        result["_fasta_file"] = self.fasta_file
        result["_genbank_file"] = self.genbank_file
        result["_protein_file"] = self.protein_file
        result["_gene_file"] = self.gene_file

        # Add resolved accessions if available
        if self.record.resolved_accessions:
            result["_resolved_accessions"] = ",".join(self.record.resolved_accessions)

        return result
