import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_inputfile(input_path: Path) -> list[tuple[str, str]]:
    """
    Parse custom input format.

    Format: First line is '<count> <length>', subsequent lines are '<name> <sequence>'.

    Returns:
        List of (organism_name, sequence) tuples.
    """
    sequences = []
    with open(input_path, encoding="utf-8") as f:
        lines = f.read().strip().splitlines()

    MIN_PARTS = 2
    # Skip header line (count + length)
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) >= MIN_PARTS:
            name = parts[0].strip()
            sequence = parts[-1].strip()
            sequences.append((name, sequence))

    return sequences


def generate_genome_accession(organism_name: str) -> str:
    """
    Convert organism name like 'P.abyssi' to accession ID like 'NC_P_abyssi'.

    Handles variations like 'Halobacterium.sp', 'M.jannaschii', etc.
    """
    # Replace dots with underscores and sanitize
    sanitized = organism_name.replace(".", "_")
    return f"NC_{sanitized}"


def prepare_raw_genomes(
    input_dir: Path,
    raw_genomes_dir: Path,
    category: str = "unclassified",
) -> dict[str, str]:
    """
    Parse inputfile and create .faa files in expected directory structure.

    Args:
        input_dir: Directory containing 'inputfile'
        raw_genomes_dir: Output directory for raw genome files
        category: Taxon category to use (default: 'unclassified')

    Returns:
        Mapping of {accession: original_name} for later restoration.
    """
    inputfile = input_dir / "inputfile"
    if not inputfile.exists():
        raise FileNotFoundError(f"Input file not found: {inputfile}")

    sequences = parse_inputfile(inputfile)
    logger.info(f"Parsed {len(sequences)} sequences from inputfile")

    # Create category directory
    category_dir = raw_genomes_dir / category
    category_dir.mkdir(parents=True, exist_ok=True)

    name_mapping: dict[str, str] = {}

    for organism_name, sequence in sequences:
        accession = generate_genome_accession(organism_name)
        name_mapping[accession] = organism_name

        # Create .faa file
        faa_path = category_dir / f"{accession}_proteins.faa"
        with open(faa_path, "w", encoding="utf-8") as f:
            f.write(f">{accession}_1 protein\n")
            f.write(f"{sequence}\n")

    logger.info(f"Created {len(name_mapping)} .faa files in {category_dir}")
    return name_mapping
