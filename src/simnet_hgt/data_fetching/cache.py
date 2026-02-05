"""
Cache manager for genome data files.

Checks if genome data already exists in the output directory before fetching.
"""

import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class GenomeCache:
    """Manages caching of genome data to avoid re-downloading existing files."""

    def __init__(self, output_dir: Path):
        """
        Initialize cache manager.

        Args:
            output_dir: Base output directory where genome files are stored
        """
        self.output_dir = Path(output_dir)
        self._cache: dict[str, set[str]] = {}  # {category: {accession1, accession2, ...}}
        logger.info(f"GenomeCache initialized for {output_dir}")

    def scan_directory(self, category: str) -> None:
        """
        Scan directory for existing genome files and populate cache.

        Args:
            category: Category directory to scan (e.g., 'archaea', 'bacteria')
        """
        category_dir = self.output_dir / category
        
        if not category_dir.exists():
            logger.debug(f"Category directory does not exist: {category_dir}")
            self._cache[category] = set()
            return

        # Look for .fasta files (the primary identifier)
        accessions = set()
        for fasta_file in category_dir.glob("*.fasta"):
            # Extract accession from filename (e.g., "NC_003326.1.fasta" -> "NC_003326.1")
            accession = fasta_file.stem
            accessions.add(accession)

        self._cache[category] = accessions
        logger.info(f"Scanned {category}: found {len(accessions)} cached accessions")

    def is_cached(self, accession: str, category: str) -> bool:
        """
        Check if a genome accession is already cached.

        Args:
            accession: RefSeq accession to check (e.g., "NC_003326.1")
            category: Category to check in

        Returns:
            True if the accession exists in cache, False otherwise
        """
        # Ensure category has been scanned
        if category not in self._cache:
            self.scan_directory(category)

        return accession in self._cache[category]

    def check_required_files(self, accession: str, category: str, 
                            check_proteins: bool = False, 
                            check_genes: bool = False) -> dict[str, bool]:
        """
        Check which required files exist for an accession.

        Args:
            accession: RefSeq accession to check
            category: Category to check in
            check_proteins: Whether to check for protein file
            check_genes: Whether to check for gene file

        Returns:
            Dictionary with file existence status
        """
        category_dir = self.output_dir / category
        
        files_status = {
            'fasta': (category_dir / f"{accession}.fasta").exists(),
            'genbank': (category_dir / f"{accession}.gb").exists(),
        }

        if check_proteins:
            files_status['proteins'] = (category_dir / f"{accession}_proteins.faa").exists()
        
        if check_genes:
            files_status['genes'] = (category_dir / f"{accession}_genes.json").exists()

        return files_status

    def has_all_required_files(self, accession: str, category: str,
                               check_proteins: bool = False,
                               check_genes: bool = False) -> bool:
        """
        Check if all required files exist for an accession.

        Args:
            accession: RefSeq accession to check
            category: Category to check in
            check_proteins: Whether protein file is required
            check_genes: Whether gene file is required

        Returns:
            True if all required files exist, False otherwise
        """
        files_status = self.check_required_files(accession, category, check_proteins, check_genes)
        
        # Core files must always exist
        has_core_files = files_status['fasta'] and files_status['genbank']
        has_required_proteins = not check_proteins or files_status.get('proteins', False)
        has_required_genes = not check_genes or files_status.get('genes', False)
        
        return has_core_files and has_required_proteins and has_required_genes

    def add_to_cache(self, accession: str, category: str) -> None:
        """
        Add an accession to the cache after successful download.

        Args:
            accession: RefSeq accession that was downloaded
            category: Category it belongs to
        """
        if category not in self._cache:
            self._cache[category] = set()
        
        self._cache[category].add(accession)
        logger.debug(f"Added {accession} to {category} cache")

    def get_cache_stats(self, category: str) -> dict[str, Any]:
        """
        Get cache statistics for a category.

        Args:
            category: Category to get stats for

        Returns:
            Dictionary with cache statistics
        """
        if category not in self._cache:
            self.scan_directory(category)

        return {
            'category': category,
            'cached_count': len(self._cache[category]),
            'cached_accessions': sorted(self._cache[category])
        }

    def clear_cache(self, category: str | None = None) -> None:
        """
        Clear the in-memory cache.

        Args:
            category: Specific category to clear, or None to clear all
        """
        if category:
            self._cache[category] = set()
            logger.info(f"Cleared cache for {category}")
        else:
            self._cache.clear()
            logger.info("Cleared all cache")
