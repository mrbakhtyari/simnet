import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)


def prepare_phylogeny(
    input_dir: Path,
    phylogeny_dir: Path,
    name_mapping: dict[str, str],
) -> None:
    """
    Process newick tree and create genome_to_organism_name.csv.

    Args:
        input_dir: Directory containing 'tree.newick'
        phylogeny_dir: Output directory for phylogeny files
        name_mapping: Mapping of {accession: original_name}
    """
    input_tree = input_dir / "tree.newick"
    if not input_tree.exists():
        raise FileNotFoundError(f"Tree file not found: {input_tree}")

    phylogeny_dir.mkdir(parents=True, exist_ok=True)

    # Read original newick tree
    with open(input_tree, encoding="utf-8") as f:
        newick_content = f.read()

    # Create reverse mapping (original_name -> accession)
    reverse_mapping = {v: k for k, v in name_mapping.items()}

    # Replace organism names with accessions in the tree
    # Sort by length descending to avoid partial replacements
    sorted_names = sorted(reverse_mapping.keys(), key=len, reverse=True)
    modified_tree = newick_content

    for original_name in sorted_names:
        accession = reverse_mapping[original_name]
        # Use word boundaries to avoid partial matches
        pattern = re.escape(original_name)
        modified_tree = re.sub(
            rf"(?<![A-Za-z0-9_.]){pattern}(?![A-Za-z0-9_])",
            accession,
            modified_tree,
        )

    # Write modified tree
    output_tree = phylogeny_dir / "tree.newick"
    with open(output_tree, "w", encoding="utf-8") as f:
        f.write(modified_tree)
    logger.info(f"Wrote modified tree to {output_tree}")

    # Create genome_to_organism_name.csv with unique taxids
    csv_path = phylogeny_dir / "genome_to_organism_name.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("genome_accession,organism_name,taxid\n")
        for idx, (accession, original_name) in enumerate(name_mapping.items(), start=1):
            f.write(f"{accession},{original_name},{idx}\n")
    logger.info(f"Wrote {len(name_mapping)} entries to {csv_path}")
