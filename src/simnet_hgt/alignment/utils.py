import csv
import logging
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def _convert_to_tsv(
    query_db: Path,
    target_db: Path,
    results_db: Path,
    output_tsv: Path,
    header: list[str],
) -> None:
    """Convert MMseqs2 binary results to TSV format."""
    logger.info("Converting results to TSV")

    cmd = [
        "mmseqs",
        "convertalis",
        str(query_db),
        str(target_db),
        str(results_db),
        str(output_tsv),
        "--format-output",
        ",".join(header),
    ]

    try:
        subprocess.run(cmd, check=True, shell=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"MMseqs2 convertalis failed: {e.stderr}") from e


def _filter_self_hits(input_tsv: Path, output_tsv: Path, header: list[str]) -> int:
    """Filter out self-hits where query == target."""
    logger.info("Filtering self-hits")
    edge_count = 0

    with (
        open(input_tsv, encoding="utf-8") as f_in,
        open(output_tsv, "w", encoding="utf-8", newline="") as f_out,
    ):
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(header)

        for row in reader:
            if row and row[0] != row[1]:  # Skip empty rows and self-hits
                writer.writerow(row)
                edge_count += 1

    return edge_count


def _clean_mutual_hits(input_tsv: Path, output_tsv: Path, header: list[str]) -> int:
    """
    De-duplicate reciprocal hits, keeping only the highest pident for each pair.

    For mutual hits (A->B and B->A), keeps only the hit with the highest 'pident'.
    Uses pandas for efficient processing of large datasets (80M+ rows).

    Args:
        input_tsv: Path to the input TSV file with protein hits
        output_tsv: Path to save the cleaned TSV file
        header: Column names for the TSV file

    Returns:
        int: Number of unique edges in the cleaned file
    """

    logger.info("Cleaning mutual hits, keeping highest pident for each pair")

    df = pd.read_csv(input_tsv, sep="\t", names=header, skiprows=1)

    # Create canonical pairs (alphabetically sorted)
    df["node1"] = df[["query", "target"]].min(axis=1)
    df["node2"] = df[["query", "target"]].max(axis=1)

    # Group by canonical pairs and keep row with highest pident
    # Using idxmax() to get the index of the max pident, then loc to select that row
    df_cleaned = df.loc[df.groupby(["node1", "node2"])["pident"].idxmax()]

    # Drop the temporary columns
    df_cleaned = df_cleaned[header]

    # Write to TSV
    df_cleaned.to_csv(output_tsv, sep="\t", index=False, header=True)

    return len(df_cleaned)


def _cleanup_temp_files(tmp_dir: Path) -> None:
    """Clean up temporary MMseqs2 files."""
    try:
        for pattern in ["search_results*", "all_hits_unfiltered.tsv", "mmseqs_db*"]:
            for f in tmp_dir.glob(pattern):
                f.unlink()

        if tmp_dir.exists() and not any(tmp_dir.iterdir()):
            tmp_dir.rmdir()
    except OSError as e:
        logger.warning("Could not clean up temporary files: %s", e)
