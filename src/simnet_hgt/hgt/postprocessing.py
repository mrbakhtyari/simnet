import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def restore_original_names(
    hgt_output_dir: Path,
    name_mapping: dict[str, str],
) -> None:
    """
    Post-process HGT results to restore original organism names.

    Creates hgt_scores_original_names.tsv with just the original organism names
    (e.g., 'P.abyssi') instead of full protein IDs.
    """
    hgt_scores_path = hgt_output_dir / "hgt_scores.tsv"
    if not hgt_scores_path.exists():
        logger.warning(f"HGT scores file not found: {hgt_scores_path}")
        return

    output_path = hgt_output_dir / "hgt_scores_original_names.tsv"

    # Create reverse mapping (accession -> original_name)
    # Used to extract original name from protein ID
    accession_to_name = name_mapping

    MIN_COLUMNS = 2

    def extract_original_name(protein_id: str) -> str:
        """Extract original organism name from protein ID."""
        # Protein ID format: category_accession_accession_1_product_1
        # e.g., unclassified_NC_P_abyssi_NC_P_abyssi_1_protein_1
        for accession, original_name in accession_to_name.items():
            if accession in protein_id:
                return original_name
        return protein_id  # Fallback to original if not found

    with open(hgt_scores_path, encoding="utf-8") as f:
        lines = f.readlines()

    with open(output_path, "w", encoding="utf-8") as f:
        # Write header
        f.write(lines[0])

        # Process data lines
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= MIN_COLUMNS:
                # Replace protein_A and protein_B with original names
                parts[0] = extract_original_name(parts[0])
                parts[1] = extract_original_name(parts[1])
            f.write("\t".join(parts) + "\n")

    logger.info(f"Wrote results with original names to {output_path}")

    # Also process high_confidence_hgt.tsv if it exists
    high_conf_path = hgt_output_dir / "high_confidence_hgt.tsv"
    if high_conf_path.exists():
        high_conf_output = hgt_output_dir / "high_confidence_hgt_original_names.tsv"

        with open(high_conf_path, encoding="utf-8") as f:
            lines = f.readlines()

        with open(high_conf_output, "w", encoding="utf-8") as f:
            # Write header
            f.write(lines[0])

            # Process data lines
            for line in lines[1:]:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= MIN_COLUMNS:
                    # Replace protein_A and protein_B with original names
                    parts[0] = extract_original_name(parts[0])
                    parts[1] = extract_original_name(parts[1])
                f.write("\t".join(parts) + "\n")

        logger.info(f"Wrote high confidence results with original names to {high_conf_output}")
