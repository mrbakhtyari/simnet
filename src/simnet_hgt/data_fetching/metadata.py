"""
GenBank metadata extractor.
"""

import logging
from typing import Any

from Bio.SeqUtils import gc_fraction

logger = logging.getLogger(__name__)


class GenBankMetadata:
    """Extracts metadata from BioPython GenBank records to match reference Excel formats."""

    @staticmethod
    def extract_all(gb_record: Any, accession: str, category: str) -> dict[str, Any]:
        """
        Extract metadata from GenBank record matching reference Excel structure.

        Args:
            gb_record: BioPython GenBank record
            accession: Accession number
            category: Dataset category (archaea, bacteria, plasmid, virus)

        Returns:
            Dictionary with category-specific metadata fields
        """
        # Category-specific extraction
        if category == "plasmid":
            return GenBankMetadata._extract_plasmid_metadata(gb_record, accession)
        elif category == "virus":
            return GenBankMetadata._extract_virus_metadata(gb_record, accession)
        elif category in ("archaea", "bacteria"):
            return GenBankMetadata._extract_archaea_bacteria_metadata(
                gb_record, accession, category
            )
        else:
            # Fallback
            return GenBankMetadata._extract_basic_metadata(gb_record, accession)

    @staticmethod
    def _extract_plasmid_metadata(gb_record: Any, accession: str) -> dict[str, Any]:
        """
        Extract plasmid metadata matching CorelS3.xlsx format.

        Required fields: UniqID, Organism/Name, TaxID, RefSeq, GI, GenBank, HostName
        """
        metadata = {
            "RefSeq": accession,
            "Organism/Name": gb_record.annotations.get("organism", ""),
            "TaxID": None,
            "GI": None,
            "GenBank": None,
            "HostName": None,
        }

        # Extract TaxID
        for feature in gb_record.features:
            if feature.type == "source" and "db_xref" in feature.qualifiers:
                for xref in feature.qualifiers["db_xref"]:
                    if xref.startswith("taxon:"):
                        metadata["TaxID"] = xref.replace("taxon:", "")
                        break
                break

        # Extract GI from annotations
        if hasattr(gb_record, "annotations"):
            gi = gb_record.annotations.get("gi")
            if gi:
                metadata["GI"] = gi

        # Extract GenBank accession from LOCUS or accession
        if hasattr(gb_record, "id"):
            metadata["GenBank"] = gb_record.id

        # Extract HostName from source feature
        for feature in gb_record.features:
            if feature.type == "source":
                if "host" in feature.qualifiers:
                    metadata["HostName"] = feature.qualifiers["host"][0]
                break

        return metadata

    @staticmethod
    def _extract_virus_metadata(gb_record: Any, accession: str) -> dict[str, Any]:
        """
        Extract virus metadata matching CorelS4.xlsx format.

        Required fields: UniqID, Organism/Name, TaxID, BioProject Accession, BioProject ID,
        Group, SubGroup, GC%, Genes, Proteins, Release Date, Modify Date, Status, Size (Kb), Host
        """
        metadata = {
            "Organism/Name": gb_record.annotations.get("organism", ""),
            "TaxID": None,
            "BioProject Accession": None,
            "BioProject ID": None,
            "Group": None,
            "SubGroup": None,
            "GC%": round(gc_fraction(gb_record.seq) * 100, 1),
            "Genes": None,  # Will be counted by parser
            "Proteins": None,  # Will be counted by parser
            "Release Date": gb_record.annotations.get("date", ""),
            "Modify Date": None,
            "Status": "Complete",
            "Size (Kb)": round(len(gb_record.seq) / 1000, 2),
            "Host": None,
        }

        # Extract identifiers and taxonomy
        GenBankMetadata._extract_taxid(gb_record, metadata)
        GenBankMetadata._extract_bioproject(gb_record, metadata)
        GenBankMetadata._extract_taxonomy_simple(gb_record, metadata)
        GenBankMetadata._extract_host(gb_record, metadata)

        return metadata

    @staticmethod
    def _extract_archaea_bacteria_metadata(
        gb_record: Any, accession: str, category: str
    ) -> dict[str, Any]:
        """
        Extract archaea/bacteria metadata matching CorelS1/CorelS2.xlsx format.

        Required fields: UniqID, Organism/Name, TaxID, BioProject Accession, BioProject ID,
        Group, SubGroup, Size (Mb), GC%, Chromosomes/RefSeq, Chromosomes/INSDC, Plasmids/RefSeq,
        Plasmids/INSDC, WGS, Scaffolds, Genes, Proteins, Release Date, Modify Date, Status, Center,
        BioSample Accession, Assembly Accession, and many phenotypic fields

        Note: Chromosomes/RefSeq and Plasmids columns will be set later from original data
        """
        metadata = {
            "Organism/Name": gb_record.annotations.get("organism", ""),
            "TaxID": None,
            "BioProject Accession": None,
            "BioProject ID": None,
            "Group": None,
            "SubGroup": None,
            "Size (Mb)": round(len(gb_record.seq) / 1_000_000, 4),
            "GC%": round(gc_fraction(gb_record.seq) * 100, 1),
            # These will be populated from original Excel data, not GenBank
            "Chromosomes/RefSeq": None,
            "Chromosomes/INSDC": None,
            "Plasmids/RefSeq": "-",
            "Plasmids/INSDC": "-",
            "WGS": "-",
            "Scaffolds": 1,
            "Genes": None,  # Will be counted by parser
            "Proteins": None,  # Will be counted by parser
            "Release Date": gb_record.annotations.get("date", ""),
            "Modify Date": None,
            "Status": "Complete",
            "Center": None,
            "BioSample Accession": "-",
            "Assembly Accession": "-",
            # Phenotypic fields - all default to None or "-"
            "Body Products": "-",
            "Biotic Relationship": None,
            "Viral Subgroup": "-",
            "Temperature Optimum": None,
            "Symbiotic Physical Interaction": "-",
            "Symbiont Name": "-",
            "Cell Length": "-",
            "Motility": None,
            "Commercial Strain Comments": "-",
            "Gram Stain": None,
            "Ph": None,
            "Oxygen Requirement": None,
            "Metabolism": None,
            "Diseases": None,
            "Commercial Strain": "-",
            "Phenotypes": "-",
            "Cell Diameter": "-",
            "Cell Arrangements": "-",
            "Temperature Range": None,
            "Sporulation": None,
            "Energy Sources": None,
            "Salinity": "-",
            "Habitats": None,
            "Symbiont Taxon ID": "-",
            "Color": "-",
            "Symbiotic Relationship": "-",
            "Viral Group": "-",
            "Cell Shape": None,
            "Kingdom": "archae" if category == "archaea" else "bacteria",  # Match reference format
            "Phylum": None,
        }

        # Extract using helper methods
        GenBankMetadata._extract_taxid(gb_record, metadata)
        GenBankMetadata._extract_bioproject_biosample(gb_record, metadata)
        # Store INSDC accession temporarily for potential use
        metadata["_current_insdc"] = gb_record.id if hasattr(gb_record, "id") else None
        GenBankMetadata._extract_taxonomy_full(gb_record, metadata, category)
        GenBankMetadata._extract_phenotypic_fields(gb_record, metadata)

        return metadata

    @staticmethod
    def _extract_basic_metadata(gb_record: Any, accession: str) -> dict[str, Any]:
        """Fallback: extract basic metadata for unknown categories."""
        metadata = {
            "Accession": accession,
            "Organism/Name": gb_record.annotations.get("organism", ""),
            "TaxID": None,
            "Size (Mb)": round(len(gb_record.seq) / 1_000_000, 2),
            "GC%": round(gc_fraction(gb_record.seq) * 100, 2),
            "Release Date": gb_record.annotations.get("date", ""),
        }

        # Extract TaxID
        for feature in gb_record.features:
            if feature.type == "source" and "db_xref" in feature.qualifiers:
                for xref in feature.qualifiers["db_xref"]:
                    if xref.startswith("taxon:"):
                        metadata["TaxID"] = xref.replace("taxon:", "")
                        break
                break

        return metadata

    # Helper methods used by extract functions

    @staticmethod
    def _extract_taxid(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract TaxID from source feature."""
        for feature in gb_record.features:
            if feature.type == "source" and "db_xref" in feature.qualifiers:
                for xref in feature.qualifiers["db_xref"]:
                    if xref.startswith("taxon:"):
                        taxid_str = xref.replace("taxon:", "")
                        # Convert to integer to match reference Excel format
                        try:
                            metadata["TaxID"] = int(taxid_str)
                        except ValueError:
                            metadata["TaxID"] = taxid_str
                        return

    @staticmethod
    def _extract_bioproject(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract BioProject from annotations."""
        for xref in gb_record.annotations.get("db_xref", []):
            if ":" in xref:
                db, value = xref.split(":", 1)
                if db == "BioProject":
                    metadata["BioProject Accession"] = value
                    if value.startswith("PRJNA"):
                        # Extract numeric ID as integer
                        try:
                            metadata["BioProject ID"] = int(value.replace("PRJNA", ""))
                        except ValueError:
                            metadata["BioProject ID"] = value.replace("PRJNA", "")

    @staticmethod
    def _extract_taxonomy_simple(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract taxonomy for Group/SubGroup (virus)."""
        taxonomy = gb_record.annotations.get("taxonomy", [])
        if len(taxonomy) > 0:
            metadata["Group"] = taxonomy[0]
        if len(taxonomy) > 1:
            metadata["SubGroup"] = taxonomy[1]

    @staticmethod
    def _extract_host(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract Host from source feature (for virus)."""
        for feature in gb_record.features:
            if feature.type == "source":
                if "host" in feature.qualifiers:
                    metadata["Host"] = feature.qualifiers["host"][0]
                return

    @staticmethod
    def _extract_bioproject_biosample(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract BioProject and BioSample from annotations."""
        for xref in gb_record.annotations.get("db_xref", []):
            if ":" not in xref:
                continue
            db, value = xref.split(":", 1)
            if db == "BioProject":
                metadata["BioProject Accession"] = value
                if value.startswith("PRJNA"):
                    # Extract numeric ID as integer
                    try:
                        metadata["BioProject ID"] = int(value.replace("PRJNA", ""))
                    except ValueError:
                        metadata["BioProject ID"] = value.replace("PRJNA", "")
            elif db == "BioSample":
                metadata["BioSample Accession"] = value

    @staticmethod
    def _extract_insdc_accession(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract INSDC accession (GenBank accession)."""
        if hasattr(gb_record, "id"):
            metadata["Chromosomes/INSDC"] = gb_record.id

    @staticmethod
    def _extract_taxonomy_full(gb_record: Any, metadata: dict[str, Any], category: str) -> None:
        """Extract taxonomy for archaea/bacteria."""
        taxonomy = gb_record.annotations.get("taxonomy", [])
        if len(taxonomy) > 0:
            metadata["Group"] = taxonomy[0]
            # Extract Phylum - typically at index 1 or 2
            # Standard taxonomy: [Kingdom, Phylum, Class, Order, Family, Genus]
            # But NCBI sometimes includes intermediate ranks
            # Phylum names often end in '-ota' or '-eota' for bacteria/archaea
            phylum = None
            for i, taxon in enumerate(taxonomy):
                if i > 0 and (taxon.endswith("ota") or taxon.endswith("eota") or i == 1):
                    phylum = taxon
                    break
            if not phylum and len(taxonomy) > 1:
                phylum = taxonomy[1]  # Fallback to second element
            metadata["Phylum"] = phylum

        if len(taxonomy) > 1:
            metadata["SubGroup"] = taxonomy[1]

    @staticmethod
    def _extract_phenotypic_fields(gb_record: Any, metadata: dict[str, Any]) -> None:
        """Extract phenotypic fields from source qualifiers."""
        for feature in gb_record.features:
            if feature.type == "source":
                qualifiers = feature.qualifiers

                # Extract available phenotypic fields
                if "isolation_source" in qualifiers:
                    metadata["Habitats"] = qualifiers["isolation_source"][0]

                if "host" in qualifiers:
                    metadata["Symbiont Name"] = qualifiers["host"][0]

                return
