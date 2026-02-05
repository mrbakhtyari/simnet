"""
GenBank file parser.
"""

import logging
import warnings
from io import StringIO
from typing import Any

from Bio import BiopythonWarning, SeqIO

logger = logging.getLogger(__name__)


class GenBankRecordParser:
    """Parses GenBank files and extracts features."""

    @staticmethod
    def parse_record(genbank_data: str) -> Any:
        """Parse GenBank text to BioPython record."""
        if not genbank_data or not genbank_data.strip():
            raise ValueError("Empty GenBank data received")

        handle = StringIO(genbank_data)
        try:
            return SeqIO.read(handle, "genbank")
        except ValueError as e:
            if "No records found" in str(e):
                raise ValueError(f"Invalid or empty GenBank format: {e!s}") from e
            raise

    @staticmethod
    def extract_proteins(gb_record: Any) -> list[str]:
        """
        Extract protein FASTA sequences from CDS features.

        Returns list of FASTA entries (one per protein).
        """
        proteins = []

        for feature in gb_record.features:
            if feature.type != "CDS":
                continue

            # Try to get pre-translated protein
            if "translation" in feature.qualifiers:
                protein_seq = feature.qualifiers["translation"][0]
            else:
                # Fallback: translate from DNA
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=BiopythonWarning)
                        protein_seq = feature.extract(gb_record.seq).translate(to_stop=True)
                except Exception:
                    logger.debug(f"Could not translate CDS at {feature.location}")
                    continue

            # Build FASTA entry
            protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0]
            product = feature.qualifiers.get("product", ["hypothetical protein"])[0]
            fasta_entry = f">{protein_id} {product}\n{protein_seq}\n"
            proteins.append(fasta_entry)

        return proteins

    @staticmethod
    def extract_genes(gb_record: Any) -> list[dict[str, Any]]:
        """
        Extract gene information from GenBank record.

        Prefers 'gene' features, falls back to 'CDS' features.
        """
        genes = []

        # Try gene features first
        gene_features = [f for f in gb_record.features if f.type == "gene"]

        if gene_features:
            for feature in gene_features:
                genes.append(
                    {
                        "type": "gene",
                        "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                        "gene": feature.qualifiers.get("gene", [""])[0],
                        "location": str(feature.location),
                    }
                )
        else:
            # Fallback to CDS (common in many genomes)
            for feature in gb_record.features:
                if feature.type == "CDS":
                    genes.append(
                        {
                            "type": "CDS",
                            "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                            "gene": feature.qualifiers.get("gene", [""])[0],
                            "protein_id": feature.qualifiers.get("protein_id", [""])[0],
                            "product": feature.qualifiers.get("product", [""])[0],
                            "location": str(feature.location),
                        }
                    )

        return genes
