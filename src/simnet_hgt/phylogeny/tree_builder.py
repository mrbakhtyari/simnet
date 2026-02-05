"""
Species tree building service.

Build phylogenetic species trees from genome taxonomy data.
"""

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez, Phylo
from Bio.Phylo.BaseTree import Clade, Tree

logger = logging.getLogger(__name__)


@dataclass
class TaxonomyInfo:
    """Taxonomy information for a single taxon."""

    taxid: int
    name: str
    lineage: str
    rank: str
    lineage_list: list[str]


def build_tree_from_lineages(
    taxonomy_data: dict,
    accession_to_taxid: dict[str, int] | None = None,
) -> Tree:
    """
    Build a phylogenetic tree from taxonomy lineages.
    Uses common ancestry to create the tree structure.

    Args:
        taxonomy_data: Dictionary mapping taxid to taxonomy info
                      Expected structure: {taxid: {"name": str, "lineage_list": list[str]}}
        accession_to_taxid: Optional mapping from genome accession to taxid.
                           If provided, leaf names will be genome accessions.
                           If None, leaf names will be "organism_name_taxid".

    Returns:
        Biopython Tree object
    """
    logger.info("Building tree from taxonomy lineages...")

    # Create root node
    root = Clade(name="Root", branch_length=0)
    tree = Tree(root=root, rooted=True)

    # Create reverse mapping: taxid -> list of accessions (if provided)
    taxid_to_accessions: dict[int, list[str]] = {}
    if accession_to_taxid:
        for acc, tid in accession_to_taxid.items():
            taxid_to_accessions.setdefault(tid, []).append(acc)

    # For each taxid, add it to the tree following its lineage
    for taxid, info in taxonomy_data.items():
        # Handle both dict and dataclass
        if hasattr(info, "lineage_list"):
            lineage_list = info.lineage_list
            scientific_name = info.name
        else:
            lineage_list = info["lineage_list"]
            scientific_name = info["name"]

        # Start from root and traverse/create nodes following the lineage
        current_node = root

        for ancestor in lineage_list:
            # Check if this ancestor already exists as a child
            found = False
            for child in current_node.clades:
                if child.name == ancestor:
                    current_node = child
                    found = True
                    break

            # If not found, create a new internal node
            if not found:
                new_node = Clade(name=ancestor, branch_length=1)
                current_node.clades.append(new_node)
                current_node = new_node

        # Add leaf node(s) for this species
        # If we have accessions for this taxid, create one leaf per accession
        if taxid in taxid_to_accessions:
            for accession in taxid_to_accessions[taxid]:
                leaf = Clade(name=accession, branch_length=1)
                current_node.clades.append(leaf)
        else:
            # Fallback: use organism_name_taxid format
            leaf = Clade(name=f"{scientific_name}_{taxid}", branch_length=1)
            current_node.clades.append(leaf)

    logger.info("Tree construction complete.")
    return tree


def fetch_taxids_from_accessions(
    accessions: list[str],
    email: str | None = None,
    api_key: str | None = None,
) -> dict[str, int]:
    """
    Fetch TaxIDs from genome accessions using NCBI esummary.
    This is a lightweight call that only gets TaxIDs, not full taxonomy data.

    Args:
        accessions: List of genome accessions (e.g., ['NC_002556.1', ...])
        email: NCBI email (required for API access)
        api_key: NCBI API key (optional, increases rate limit)

    Returns:
        Dictionary mapping accession to taxid
    """
    if email:
        Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    logger.info(f"Fetching TaxIDs for {len(accessions)} genome accessions...")

    accession_to_taxid: dict[str, int] = {}
    batch_size = 200  # esummary is lightweight, can use larger batches

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        ids_str = ",".join(batch)

        try:
            handle = Entrez.esummary(db="nucleotide", id=ids_str, retmode="xml")
            summaries = Entrez.read(handle)
            handle.close()

            # Handle the response - summaries can be a list or dict
            if isinstance(summaries, list):
                summary_list = summaries
            else:
                # Sometimes it's a dict with DocSum entries
                summary_list = [summaries] if summaries else []

            for summary in summary_list:
                # Type guard - summary should be a dict-like object
                if not isinstance(summary, dict):
                    continue

                accession = summary.get("AccessionVersion", summary.get("Caption", ""))
                taxid_str = summary.get("TaxId", "0")

                if accession and taxid_str:
                    taxid = int(taxid_str)
                    if taxid > 0:
                        accession_to_taxid[accession] = taxid

        except Exception as e:
            logger.warning(f"Error fetching TaxIDs for batch starting at {i}: {e}")
            continue

    logger.info(f"Mapped {len(accession_to_taxid)} accessions to TaxIDs")
    return accession_to_taxid


def fetch_taxonomy_info(
    taxid_list: list[int],
    email: str | None = None,
    api_key: str | None = None,
) -> dict[int, TaxonomyInfo]:
    """
    Fetch taxonomy information from NCBI for a list of TaxIDs.

    Args:
        taxid_list: List of NCBI taxonomy IDs
        email: NCBI email (required for API access)
        api_key: NCBI API key (optional, increases rate limit)

    Returns:
        Dictionary mapping taxid to TaxonomyInfo
    """
    if email:
        Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    logger.info(f"Fetching taxonomy info for {len(taxid_list)} TaxIDs from NCBI...")

    taxonomy_data: dict[int, TaxonomyInfo] = {}
    batch_size = 100  # NCBI allows batching

    for i in range(0, len(taxid_list), batch_size):
        batch = taxid_list[i : i + batch_size]
        ids_str = ",".join(str(tid) for tid in batch)

        try:
            handle = Entrez.efetch(db="taxonomy", id=ids_str, retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            for record in records:  # type: ignore
                taxid = int(record["TaxId"])
                scientific_name = record.get("ScientificName", f"Unknown_{taxid}")
                lineage = record.get("Lineage", "")
                rank = record.get("Rank", "")
                lineage_list = [item.strip() for item in lineage.split(";") if item.strip()]

                taxonomy_data[taxid] = TaxonomyInfo(
                    taxid=taxid,
                    name=scientific_name,
                    lineage=lineage,
                    rank=rank,
                    lineage_list=lineage_list,
                )

        except Exception as e:
            logger.warning(f"Error fetching batch starting at index {i}: {e}")
            continue

    logger.info(f"Successfully fetched info for {len(taxonomy_data)} TaxIDs.")
    return taxonomy_data


def build_species_tree(
    genome_accessions: list[str],
    output_path: Path,
    prefix: str = "",
    email: str | None = None,
    api_key: str | None = None,
) -> Path:
    """
    Build a species tree from a list of genome accessions.

    Args:
        genome_accessions: List of genome accessions (e.g., ['NC_002556.1', ...])
        output_path: Directory to save output files.
        prefix: Prefix for output file names.
        email: NCBI email (required for API access).
        api_key: NCBI API key (optional, increases rate limit).

    Returns:
        Path to the generated Newick tree file, or None on failure.
    """

    output_path.mkdir(parents=True, exist_ok=True)
    output_newick = output_path / f"{prefix}_species_tree.nwk"

    # Step 1: Lightweight call to get TaxIDs from accessions
    accession_to_taxid = fetch_taxids_from_accessions(
        genome_accessions, email=email, api_key=api_key
    )

    # Step 2: Single taxonomy API call for all unique TaxIDs
    unique_taxids = list(set(accession_to_taxid.values()))
    taxonomy_data = fetch_taxonomy_info(unique_taxids, email=email, api_key=api_key)

    # Build the phylogenetic tree with genome accessions as leaf names
    tree = build_tree_from_lineages(taxonomy_data, accession_to_taxid)

    # Save the tree to a Newick file
    logger.info(f"Saving tree to {output_newick}...")
    try:
        Phylo.write(tree, output_newick, "newick")  # type: ignore
        logger.info("Tree saved successfully.")
    except Exception as e:
        logger.error(f"Error writing tree file: {e}")

    # Save genome accession to organism name mapping for reference
    # Note: Tree now uses genome accessions as leaf names, so this is optional metadata
    genome_to_organism_file = output_newick.parent / "genome_to_organism_name.csv"
    logger.info(f"Saving genome-to-organism mapping to {genome_to_organism_file}...")
    try:
        with open(genome_to_organism_file, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["genome_accession", "organism_name", "taxid"])
            for accession, taxid in accession_to_taxid.items():
                if taxid in taxonomy_data:
                    organism_name = taxonomy_data[taxid].name
                    writer.writerow([accession, organism_name, taxid])
        logger.info("Genome-to-organism mapping saved successfully.")
    except Exception as e:
        logger.warning(f"Could not save genome-to-organism mapping: {e}")

    # Calculate statistics
    leaf_count = sum(1 for _ in tree.get_terminals())
    total_count = tree.count_terminals() + len(list(tree.get_nonterminals()))

    logger.info("--- Species Tree Complete ---")
    logger.info(f"Tree saved in Newick format to: {output_newick}")
    logger.info(f"Number of leaves (genome accessions): {leaf_count}")
    logger.info(f"Total nodes (including internal): {total_count}")

    return output_newick
