"""
Genome fetcher service.
"""

import json
import logging
from pathlib import Path

from tqdm import tqdm

from .cache import GenomeCache
from .client import EntrezClient
from .config import DatasetConfig
from .genome import FetchResult, OrganismRecord
from .metadata import GenBankMetadata
from .parser import GenBankRecordParser
from .storage import GenomeStorage

logger = logging.getLogger(__name__)


class GenomeFetcher:
    """Orchestrates genome data fetching from NCBI."""

    def __init__(
        self,
        output_dir: Path,
        email: str | None = None,
        api_key: str | None = None,
        delay: float = 0.4,
        use_cache: bool = True,
    ):
        """
        Initialize fetcher.

        Args:
            output_dir: Base output directory
            client: NCBI client (for dependency injection / testing)
            email: NCBI email (required if client not provided)
            api_key: NCBI API key (optional, increases rate limit)
            delay: Delay between requests
            use_cache: Whether to use cache to skip already downloaded genomes
        """
        self.output_dir = Path(output_dir)
        self.use_cache = use_cache

        # Dependency injection: use provided client or create default
        self.client = EntrezClient(email=email, api_key=api_key, delay=delay)
        self.parser = GenBankRecordParser()
        self.metadata_extractor = GenBankMetadata()
        self.file_writer = GenomeStorage(output_dir)

        # Initialize cache manager
        self.cache = GenomeCache(output_dir)

        # Track fetched accessions to avoid duplicates
        self.fetched_accessions: dict[str, str] = {}  # {refseq_accession: uniq_id}

        logger.info(
            f"NCBIFetcher initialized: {output_dir} (cache={'enabled' if use_cache else 'disabled'})"
        )

    def fetch_dataset(
        self,
        config: DatasetConfig,
        max_records: int | None = None,
    ) -> None:
        """Fetch all records in a dataset."""
        logger.info(f"Fetching {config.name} dataset")

        # Scan cache for existing files if cache is enabled
        if self.use_cache:
            self.cache.scan_directory(config.name)
            cache_stats = self.cache.get_cache_stats(config.name)
            logger.info(f"Cache contains {cache_stats['cached_count']} existing accessions")

        records = config.load_records()
        if max_records:
            records = records[:max_records]

        if not records:
            logger.warning(f"No records found for {config.name}")

        results = []

        with tqdm(records, desc=f"{config.name}", unit="organism") as pbar:
            for record in pbar:
                display = record.organism or record.get_primary_accession()
                pbar.set_postfix_str(display[:40])

                result = self._fetch_organism(record, config)
                results.append(result)

        # Save metadata summary
        self.file_writer.save_metadata(results, config.name)

        success = sum(1 for r in results if r.is_success())
        cached = sum(1 for r in results if r.status == "cached")
        logger.info(
            f"Completed {config.name}: {success}/{len(records)} succeeded ({cached} from cache)"
        )

    def _fetch_organism(self, record: OrganismRecord, config: DatasetConfig) -> FetchResult:
        """Fetch data for one organism."""
        try:
            # Resolve BioProject if needed
            accessions = self._resolve_accessions(record)

            # Check cache BEFORE fetching anything from NCBI
            if self.use_cache:
                # For combined chromosomes, check cache with first accession
                if config.combine_chromosomes and len(accessions) > 1:
                    primary_accession = accessions[0]
                    if self.cache.has_all_required_files(
                        primary_accession,
                        config.name,
                        check_proteins=config.fetch_proteins,
                        check_genes=config.fetch_genes,
                    ):
                        logger.info(f"Using cached combined data for {primary_accession}")
                        return self._load_from_cache(primary_accession, record, config)
                else:
                    # For single accession, check cache
                    accession = accessions[0]
                    if not accession.startswith("PRJ") and self.cache.has_all_required_files(
                        accession,
                        config.name,
                        check_proteins=config.fetch_proteins,
                        check_genes=config.fetch_genes,
                    ):
                        logger.info(f"Using cached data for {accession}")
                        return self._load_from_cache(accession, record, config)

            # Fetch sequences (only if not cached)
            if config.combine_chromosomes and len(accessions) > 1:
                result = self._fetch_organism_combined_chromosomes(record, accessions, config)
            else:
                result = self._fetch_organism_single_accession(record, accessions[0], config)

            # Validate if requested
            if config.validate_counts:
                self._validate_counts(result, record)

            return result

        except Exception as e:
            logger.error(f"Failed {record.get_primary_accession()}: {e}")
            return FetchResult(record=record, status="failed", error=str(e))

    def _resolve_accessions(self, record: OrganismRecord) -> list[str]:
        """Resolve BioProject to accessions if needed."""
        primary = record.get_primary_accession()

        if primary.startswith("PRJ"):
            bpid = record.additional_fields.get("BioProject ID", primary)

            logger.debug(
                f"Resolving BioProject {primary} for {record.organism} using BioProject ID {bpid}"
            )

            if bpid:
                project_term = f"({primary}[BioProject] OR {bpid}[BioProject ID])"
            else:
                project_term = f"{primary}[BioProject]"

            # BioProject - try multiple search strategies
            search_terms = [
                f"{project_term} AND (complete genome[filter] OR complete sequence[filter]) AND refseq[filter]",
                f"{project_term} AND (complete genome[filter] OR complete sequence[filter])",
                f"{project_term}",
            ]

            for term in search_terms:
                try:
                    ids = self.client.search(term, retmax=1)
                    if ids:
                        record.resolved_accessions = ids
                        logger.debug(f"Resolved {primary} to {ids[0]} using: {term}")
                        return ids
                except Exception as e:
                    logger.debug(f"Search failed with '{term}': {e}")
                    continue

            # If BioProject search failed, try organism name (common for old virus projects)
            if record.organism:
                logger.debug(f"BioProject search failed, trying organism name: {record.organism}")
                organism_terms = [
                    f'"{record.organism}"[Organism] AND (complete genome[title] OR complete sequence[title]) AND refseq[filter]',
                    f'"{record.organism}"[Organism] AND (complete genome[title] OR complete sequence[title])',
                    f'"{record.organism}"[Organism] AND refseq[filter]',
                    f'"{record.organism}"[Organism]',
                ]

                for term in organism_terms:
                    try:
                        ids = self.client.search(term, retmax=1)
                        if ids:
                            record.resolved_accessions = ids
                            logger.debug(
                                f"Resolved {primary} via organism name to {ids[0]} using: {term}"
                            )
                            return ids
                    except Exception as e:
                        logger.debug(f"Organism search failed with '{term}': {e}")
                        continue

            raise ValueError(f"Could not resolve BioProject {primary}")

        return record.accessions

    def _fetch_organism_single_accession(
        self, record: OrganismRecord, accession: str, config: DatasetConfig
    ) -> FetchResult:
        """Fetch single sequence."""

        # Fetch data from NCBI
        fasta_data = self.client.fetch_sequence(accession, "fasta")
        genbank_data = self.client.fetch_sequence(accession, "gb")

        # Parse GenBank
        gb_record = self.parser.parse_record(genbank_data)

        # Get RefSeq accession (NC_XXXXXX format) from the GenBank record
        refseq_accession = gb_record.id  # This is the RefSeq accession like NC_003326.1

        # Check if the RefSeq accession is different and already cached
        if (
            self.use_cache
            and refseq_accession != accession
            and self.cache.has_all_required_files(
                refseq_accession,
                config.name,
                check_proteins=config.fetch_proteins,
                check_genes=config.fetch_genes,
            )
        ):
            logger.info(f"Using cached data for {refseq_accession} (resolved from {accession})")
            return self._load_from_cache(refseq_accession, record, config)

        # Check if we already fetched this accession from a different BioProject
        if refseq_accession in self.fetched_accessions:
            previous_uniq_id = self.fetched_accessions[refseq_accession]
            logger.info(
                f"Skipping duplicate: {refseq_accession} already fetched for {previous_uniq_id}, "
                f"now requested for {record.uniq_id}"
            )
            return FetchResult(
                record=record,
                status="duplicate",
                error=f"Already fetched as {previous_uniq_id}",
                sequence_length=0,
            )

        # Mark this accession as fetched
        self.fetched_accessions[refseq_accession] = record.uniq_id

        # Extract metadata from GenBank
        metadata = self.metadata_extractor.extract_all(gb_record, accession, config.name)
        metadata["RefSeq_Accession"] = refseq_accession  # Add RefSeq accession to metadata

        # Extract features
        proteins = self.parser.extract_proteins(gb_record) if config.fetch_proteins else []
        genes = self.parser.extract_genes(gb_record) if config.fetch_genes else []

        # Update metadata with actual counts (only for categories that have these columns)
        if config.name in ("archaea", "bacteria", "virus"):
            metadata["Genes"] = len(genes)
            metadata["Proteins"] = len(proteins)

        # For archaea/bacteria, set accession columns from original Excel data
        if config.name in ("archaea", "bacteria"):
            metadata["Chromosomes/RefSeq"] = record.additional_fields.get(
                "_original_chromosomes_refseq", accession
            )
            metadata["Chromosomes/INSDC"] = record.additional_fields.get(
                "_original_chromosomes_insdc", metadata.get("_current_insdc", "")
            )
            metadata["Plasmids/RefSeq"] = record.additional_fields.get(
                "_original_plasmids_refseq", "-"
            )
            metadata["Plasmids/INSDC"] = record.additional_fields.get(
                "_original_plasmids_insdc", "-"
            )
            # Remove temporary field
            metadata.pop("_current_insdc", None)

        # Update record with enriched metadata
        record.additional_fields.update(metadata)

        # Save files using RefSeq accession as base name
        fasta_file = self.file_writer.save_fasta(fasta_data, refseq_accession, config.name)
        gb_file = self.file_writer.save_genbank(genbank_data, refseq_accession, config.name)

        protein_file = None
        if proteins:
            protein_file = self.file_writer.save_proteins(proteins, refseq_accession, config.name)

        gene_file = None
        if genes:
            gene_file = self.file_writer.save_genes(genes, refseq_accession, config.name)

        return FetchResult(
            record=record,
            status="success",
            sequence_length=len(fasta_data),
            num_genes=len(genes),
            num_proteins=len(proteins),
            fasta_file=str(fasta_file.relative_to(self.output_dir)),
            genbank_file=str(gb_file.relative_to(self.output_dir)),
            protein_file=str(protein_file.relative_to(self.output_dir)) if protein_file else None,
            gene_file=str(gene_file.relative_to(self.output_dir)) if gene_file else None,
        )

    def _load_from_cache(
        self, accession: str, record: OrganismRecord, config: DatasetConfig
    ) -> FetchResult:
        """Load genome data from cached files."""
        category_dir = self.output_dir / config.name

        # Read cached GenBank file to extract metadata
        gb_file = category_dir / f"{accession}.gb"
        genbank_data = gb_file.read_text(encoding="utf-8")
        gb_record = self.parser.parse_record(genbank_data)

        # Extract metadata from GenBank
        metadata = self.metadata_extractor.extract_all(gb_record, accession, config.name)
        metadata["RefSeq_Accession"] = gb_record.id

        # Count genes and proteins from cached files if they exist
        num_genes = 0
        num_proteins = 0

        if config.fetch_genes:
            gene_file = category_dir / f"{accession}_genes.json"
            if gene_file.exists():
                genes = json.loads(gene_file.read_text(encoding="utf-8"))
                num_genes = len(genes)

        if config.fetch_proteins:
            protein_file = category_dir / f"{accession}_proteins.faa"
            if protein_file.exists():
                protein_data = protein_file.read_text(encoding="utf-8")
                num_proteins = protein_data.count(">")

        # Update metadata with actual counts
        if config.name in ("archaea", "bacteria", "virus"):
            metadata["Genes"] = num_genes
            metadata["Proteins"] = num_proteins

        # For archaea/bacteria, set accession columns
        if config.name in ("archaea", "bacteria"):
            metadata["Chromosomes/RefSeq"] = record.additional_fields.get(
                "_original_chromosomes_refseq", accession
            )
            metadata["Chromosomes/INSDC"] = record.additional_fields.get(
                "_original_chromosomes_insdc", metadata.get("_current_insdc", "")
            )
            metadata["Plasmids/RefSeq"] = record.additional_fields.get(
                "_original_plasmids_refseq", "-"
            )
            metadata["Plasmids/INSDC"] = record.additional_fields.get(
                "_original_plasmids_insdc", "-"
            )
            metadata.pop("_current_insdc", None)

        # Update record with metadata
        record.additional_fields.update(metadata)

        # Mark this accession as fetched (from cache)
        self.fetched_accessions[accession] = record.uniq_id

        # Get file paths
        fasta_file = category_dir / f"{accession}.fasta"
        protein_file = category_dir / f"{accession}_proteins.faa" if config.fetch_proteins else None
        gene_file = category_dir / f"{accession}_genes.json" if config.fetch_genes else None

        # Get sequence length
        sequence_length = len(fasta_file.read_text(encoding="utf-8")) if fasta_file.exists() else 0

        return FetchResult(
            record=record,
            status="cached",
            sequence_length=sequence_length,
            num_genes=num_genes,
            num_proteins=num_proteins,
            fasta_file=str(fasta_file.relative_to(self.output_dir)),
            genbank_file=str(gb_file.relative_to(self.output_dir)),
            protein_file=str(protein_file.relative_to(self.output_dir))
            if protein_file and protein_file.exists()
            else None,
            gene_file=str(gene_file.relative_to(self.output_dir))
            if gene_file and gene_file.exists()
            else None,
        )

    def _fetch_organism_combined_chromosomes(
        self, record: OrganismRecord, accessions: list[str], config: DatasetConfig
    ) -> FetchResult:
        """Fetch and combine multiple sequences (chromosomes + plasmids)."""
        logger.info(f"  Combining {len(accessions)} sequences")

        all_fasta = []
        all_genbank = []
        all_proteins = []
        all_genes = []
        all_metadata = []

        for acc in accessions:
            fasta = self.client.fetch_sequence(acc, "fasta")
            genbank = self.client.fetch_sequence(acc, "gb")
            gb_record = self.parser.parse_record(genbank)

            all_fasta.append(fasta)
            all_genbank.append(genbank)
            all_metadata.append(self.metadata_extractor.extract_all(gb_record, acc, config.name))

            if config.fetch_proteins:
                all_proteins.extend(self.parser.extract_proteins(gb_record))
            if config.fetch_genes:
                all_genes.extend(self.parser.extract_genes(gb_record))

        # Merge metadata (prioritize first accession)
        merged_metadata = all_metadata[0]

        # Sum sizes for combined genomes
        if "Size (Mb)" in merged_metadata:
            merged_metadata["Size (Mb)"] = round(
                sum(m.get("Size (Mb)", 0) for m in all_metadata), 4
            )

        # Update counts (only for categories that have these columns)
        if config.name in ("archaea", "bacteria", "virus"):
            merged_metadata["Genes"] = len(all_genes)
            merged_metadata["Proteins"] = len(all_proteins)

        # For archaea/bacteria, set accession columns from original Excel data
        if config.name in ("archaea", "bacteria"):
            merged_metadata["Chromosomes/RefSeq"] = record.additional_fields.get(
                "_original_chromosomes_refseq", ",".join(accessions)
            )
            # Collect all INSDC accessions
            insdc_list = [
                m.get("_current_insdc", "") for m in all_metadata if m.get("_current_insdc")
            ]
            merged_metadata["Chromosomes/INSDC"] = record.additional_fields.get(
                "_original_chromosomes_insdc", ",".join(insdc_list) if insdc_list else ""
            )
            merged_metadata["Plasmids/RefSeq"] = record.additional_fields.get(
                "_original_plasmids_refseq", "-"
            )
            merged_metadata["Plasmids/INSDC"] = record.additional_fields.get(
                "_original_plasmids_insdc", "-"
            )
            # Remove temporary fields
            merged_metadata.pop("_current_insdc", None)

        record.additional_fields.update(merged_metadata)

        # Save combined files
        base_name = accessions[0]
        fasta_file = self.file_writer.save_fasta("\n".join(all_fasta), base_name, config.name)
        gb_file = self.file_writer.save_genbank("\n".join(all_genbank), base_name, config.name)

        protein_file = None
        if all_proteins:
            protein_file = self.file_writer.save_proteins(all_proteins, base_name, config.name)

        gene_file = None
        if all_genes:
            gene_file = self.file_writer.save_genes(all_genes, base_name, config.name)

        return FetchResult(
            record=record,
            status="success",
            sequence_length=sum(len(f) for f in all_fasta),
            num_genes=len(all_genes),
            num_proteins=len(all_proteins),
            fasta_file=str(fasta_file.relative_to(self.output_dir)),
            genbank_file=str(gb_file.relative_to(self.output_dir)),
            protein_file=str(protein_file.relative_to(self.output_dir)) if protein_file else None,
            gene_file=str(gene_file.relative_to(self.output_dir)) if gene_file else None,
        )

    def _validate_counts(self, result: FetchResult, record: OrganismRecord) -> None:
        """Validate extracted vs. expected counts (±10%)."""
        expected_genes = record.additional_fields.get("_expected_genes")
        expected_proteins = record.additional_fields.get("_expected_proteins")

        if expected_genes and not (
            0.9 * expected_genes <= result.num_genes <= 1.1 * expected_genes
        ):
            logger.warning(f"  Gene mismatch: got {result.num_genes}, expected {expected_genes}")

        if expected_proteins and not (
            0.9 * expected_proteins <= result.num_proteins <= 1.1 * expected_proteins
        ):
            logger.warning(
                f"  Protein mismatch: got {result.num_proteins}, expected {expected_proteins}"
            )
