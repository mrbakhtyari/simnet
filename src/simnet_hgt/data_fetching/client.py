"""
NCBI Entrez API client.
"""

import logging
import time
from abc import ABC, abstractmethod
from functools import wraps
from urllib.error import HTTPError, URLError

from Bio import Entrez

logger = logging.getLogger(__name__)


def retry_on_network_error(max_retries: int = 3, backoff: float = 2.0):
    """Retry on network errors with exponential backoff."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except (OSError, HTTPError, URLError) as e:
                    last_exception = e
                    if attempt < max_retries - 1:
                        wait = backoff**attempt
                        logger.warning(f"Retry {attempt + 1}/{max_retries} in {wait}s: {e}")
                        time.sleep(wait)
            raise last_exception  # type: ignore

        return wrapper

    return decorator


class NCBIClient(ABC):
    """Abstract interface for NCBI operations."""

    @abstractmethod
    def fetch_sequence(self, accession: str, rettype: str) -> str:
        """Fetch sequence data."""
        pass

    @abstractmethod
    def search(self, term: str, retmax: int = 1) -> list[str]:
        """Search for sequences."""
        pass

    @abstractmethod
    def fetch_bioproject(self, bioproject_id: str) -> dict:
        """Fetch BioProject metadata."""
        pass

    @abstractmethod
    def fetch_biosample(self, biosample_id: str) -> dict:
        """Fetch BioSample metadata."""
        pass

    @abstractmethod
    def fetch_assembly(self, assembly_accession: str) -> dict:
        """Fetch Assembly metadata."""
        pass


class EntrezClient(NCBIClient):
    """Real NCBI Entrez client."""

    def __init__(self, email: str | None = None, api_key: str | None = None, delay: float = 0.4):
        """
        Initialize Entrez client.

        Args:
            email: Required by NCBI
            api_key: Optional, increases rate limit from 3/s to 10/s
            delay: Delay between requests (0.34s with key, 0.4s without)
        """
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self.delay = delay
        logger.info(f"EntrezClient initialized (delay={delay}s)")

    @retry_on_network_error(max_retries=3, backoff=2.0)
    def fetch_sequence(self, accession: str, rettype: str) -> str:
        """Fetch sequence from NCBI."""
        actual_rettype = "gbwithparts" if rettype == "gb" else rettype

        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype=actual_rettype, retmode="text"
        )
        data = handle.read()
        handle.close()

        time.sleep(self.delay)
        return data

    @retry_on_network_error(max_retries=3, backoff=2.0)
    def search(self, term: str, retmax: int = 1) -> list[str]:
        """Search NCBI nucleotide database."""
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()

        time.sleep(self.delay)
        return record.get("IdList", [])  # type: ignore

    @retry_on_network_error(max_retries=3, backoff=2.0)
    def fetch_bioproject(self, bioproject_id: str) -> dict:
        """Fetch BioProject metadata."""
        try:
            handle = Entrez.efetch(db="bioproject", id=bioproject_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            time.sleep(self.delay)

            if not records or "DocumentSummarySet" not in records:  # type: ignore
                return {}

            doc_sum = records["DocumentSummarySet"]["DocumentSummary"][0]  # type: ignore

            return {
                "BioProject ID": doc_sum.get("Project_Id", ""),
                "BioProject Accession": doc_sum.get("Project_Acc", bioproject_id),
                "Project Name": doc_sum.get("Project_Title", ""),
                "Project Description": doc_sum.get("Project_Description", ""),
                "Center": doc_sum.get("Registration", ""),
            }
        except Exception as e:
            logger.warning(f"Could not fetch BioProject {bioproject_id}: {e}")
            return {}

    @retry_on_network_error(max_retries=3, backoff=2.0)
    def fetch_biosample(self, biosample_id: str) -> dict:
        """Fetch BioSample metadata with phenotypic attributes."""
        try:
            handle = Entrez.efetch(db="biosample", id=biosample_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            time.sleep(self.delay)

            if not records or "BioSample" not in records:  # type: ignore
                return {}

            biosample = records["BioSample"][0]  # type: ignore

            metadata = {
                "BioSample Accession": biosample.get("@accession", biosample_id),
                "Sample Name": biosample.get("Title", ""),
            }

            # Extract attributes
            if "Attributes" in biosample:
                attrs = biosample["Attributes"]
                if isinstance(attrs, dict):
                    attrs = [attrs]

                for attr in attrs:
                    attr_name = attr.get("@attribute_name", "")
                    attr_value = str(attr.get("#text", ""))

                    # Map common attributes
                    attr_mapping = {
                        "strain": "Strain",
                        "isolation_source": "Habitats",
                        "host": "Host",
                        "collection_date": "Collection Date",
                        "geo_loc_name": "Geographic Location",
                        "temperature": "Temperature Optimum",
                        "gram_stain": "Gram Stain",
                        "cell_shape": "Cell Shape",
                        "motility": "Motility",
                        "sporulation": "Sporulation",
                        "oxygen_requirement": "Oxygen Requirement",
                        "salinity": "Salinity",
                        "pH": "Ph",
                        "phenotype": "Phenotypes",
                    }

                    if attr_name.lower() in attr_mapping:
                        metadata[attr_mapping[attr_name.lower()]] = attr_value

            return metadata

        except Exception as e:
            logger.warning(f"Could not fetch BioSample {biosample_id}: {e}")
            return {}

    @retry_on_network_error(max_retries=3, backoff=2.0)
    def fetch_assembly(self, assembly_accession: str) -> dict:
        """Fetch Assembly metadata."""
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_accession)
            records = Entrez.read(handle)
            handle.close()
            time.sleep(self.delay)

            if not records or "DocumentSummarySet" not in records:  # type: ignore
                return {}

            doc_sum = records["DocumentSummarySet"]["DocumentSummary"][0]  # type: ignore

            return {
                "Assembly Accession": doc_sum.get("AssemblyAccession", assembly_accession),
                "Assembly Name": doc_sum.get("AssemblyName", ""),
                "Status": doc_sum.get("AssemblyStatus", ""),
                "Modify Date": doc_sum.get("LastUpdateDate", ""),
                "Scaffolds": doc_sum.get("ScaffoldN50", ""),
                "WGS": doc_sum.get("WGS", ""),
            }

        except Exception as e:
            logger.warning(f"Could not fetch Assembly {assembly_accession}: {e}")
            return {}
