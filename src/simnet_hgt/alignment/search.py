"""MMseqs2 search functionality for all-vs-all protein alignment."""

import logging
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

from .utils import (
    _clean_mutual_hits,
    _cleanup_temp_files,
    _convert_to_tsv,
    _filter_self_hits,
)

logger = logging.getLogger(__name__)


@dataclass
class SearchParams:
    """Parameters for MMseqs2 search."""

    min_coverage: float = 0.8
    min_identity: float = 0.3
    coverage_mode: int = 0
    e_value: float = 1e-5


def align_sequences(
    fasta_path: Path,
    output_dir: Path,
    params: SearchParams | None = None,
    remove_self_hits: bool = True,
) -> tuple[Path, float]:
    """
    Run all-vs-all protein alignment using MMseqs2.

    Args:
        fasta_path: Path to input FASTA file
        output_dir: Directory to store output files
        params: Search parameters (min_coverage, min_identity, coverage_mode, evalue)
        remove_self_hits: If True, remove self-hits from results (default: True)

    Returns:
        Tuple containing:
        - Path to the final cleaned TSV file
        - Computation time in seconds for the search step


    Raises:
        RuntimeError: If MMseqs2 commands fail
        FileNotFoundError: If FASTA file or MMseqs2 not found
    """
    if params is None:
        params = SearchParams()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Create MMseqs2 database
    db_path = output_dir / "mmseqs_db" / "all_proteins_db"
    _create_database(fasta_path, db_path)

    # Step 2: Run search
    tmp_dir = output_dir / "tmp_search"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    search_results = tmp_dir / "search_results"
    search_time = _run_mmseqs_search(db_path, db_path, search_results, tmp_dir, params)

    # Step 3: Convert to TSV
    raw_hits_tsv = tmp_dir / "hits_raw.tsv"
    header = ["query", "target", "pident", "alnlen", "evalue", "bits"]
    _convert_to_tsv(db_path, db_path, search_results, raw_hits_tsv, header)

    # Step 4: Filter self-hits if requested
    filtered_hits_tsv = output_dir / "filtered_hits.tsv"
    if remove_self_hits:
        edge_count = _filter_self_hits(raw_hits_tsv, filtered_hits_tsv, header)
        logger.info("Filtered out self-hits, %s edges remaining", f"{edge_count:,}")
    else:
        raw_hits_tsv.rename(filtered_hits_tsv)
        with open(filtered_hits_tsv, encoding="utf-8") as f:
            edge_count = sum(1 for _ in f) - 1
        logger.info("Kept all hits including self-hits, %s edges total", f"{edge_count:,}")

    # Step 5: Deduplicate reciprocal hits (keep only highest pident for each pair)
    final_hits_tsv = output_dir / "hits_cleaned.tsv"
    unique_edge_count = _clean_mutual_hits(filtered_hits_tsv, final_hits_tsv, header)
    logger.info("Deduplicated reciprocal hits, %s unique edges saved", f"{unique_edge_count:,}")

    # Step 6: Remove intermediate file
    if filtered_hits_tsv.exists():
        filtered_hits_tsv.unlink()
        logger.debug("Removed intermediate file: %s", filtered_hits_tsv)

    # Step 7: Cleanup temporary files
    _cleanup_temp_files(tmp_dir)

    return final_hits_tsv, search_time


def _run_mmseqs_search(
    query_db: Path,
    target_db: Path,
    results_db: Path,
    tmp_dir: Path,
    params: SearchParams,
) -> float:
    """
    Execute MMseqs2 search command.

    Returns:
        Execution time in seconds.
    """
    logger.info("Running all-vs-all protein search")

    cmd = [
        "mmseqs",
        "search",
        str(query_db),
        str(target_db),
        str(results_db),
        str(tmp_dir),
        "-e",
        str(params.e_value),
        "--min-seq-id",
        str(params.min_identity),
        "-c",
        str(params.min_coverage),
        "--cov-mode",
        str(params.coverage_mode),
    ]

    try:
        start_time = time.perf_counter()
        subprocess.run(cmd, check=True, shell=True, capture_output=True, text=True)
        end_time = time.perf_counter()
        return end_time - start_time
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"MMseqs2 search failed: {e.stderr}") from e
    except FileNotFoundError as e:
        raise FileNotFoundError("MMseqs2 not found. Please install MMseqs2.") from e


def _create_database(fasta_path: Path, db_path: Path) -> None:
    """
    Create MMseqs2 database from FASTA file.

    Args:
        fasta_path: Path to input FASTA file
        db_path: Path to output MMseqs2 database

    Raises:
        RuntimeError: If MMseqs2 indexing fails
        FileNotFoundError: If FASTA file or MMseqs2 is not found
    """
    fasta_path = Path(fasta_path)
    db_path = Path(db_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    db_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        subprocess.run(
            ["mmseqs", "createdb", str(fasta_path), str(db_path)],
            check=True,
            shell=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        msg = f"MMseqs2 indexing failed: {e.stderr}"
        logger.error(msg)
        raise RuntimeError(msg) from e
    except FileNotFoundError as e:
        msg = "MMseqs2 not found. Please install MMseqs2."
        logger.error(msg)
        raise FileNotFoundError(msg) from e
