"""
Module for mapping protein hits to genomes and taxa.
"""

from pathlib import Path

import numpy as np
import pandas as pd


def create_genome_hits(
    protein_hits_path: str | Path,
    protein_metadata_path: str | Path,
    chg_membership_path: str | Path,
    output_dir: str | Path,
) -> None:
    """
    Process protein hits to filter for inter-genome hits and add taxonomic information.

    This function takes protein-to-protein BLAST hits and:
    1. Adds genome_accession for both query and target proteins
    2. Filters out hits where query and target are from the same genome
    3. Adds CHG (Co-Homologous Group) IDs for both query and target proteins
    4. Filters to keep only hits where query and target have the same CHG ID
    5. Adds taxon information for both query and target
    6. Creates a normalized taxon pair column (alphabetically sorted)
    7. Saves the result as genome_hits.tsv

    Parameters
    ----------
    protein_hits_path : str | Path
        Path to all_protein_hits.tsv containing columns:
        query, target, pident, alnlen, evalue, bits
    protein_metadata_path : str | Path
        Path to protein_metadata.csv containing columns:
        protein_id, genome_accession, taxon
    chg_membership_path : str | Path
        Path to reduced_chg_membership.tsv containing columns:
        protein_id, chg_id, original_chg_id
    output_dir : str | Path
        Directory where genome_hits.tsv will be saved

    Returns
    -------
    None
        Saves genome_hits.tsv to the output directory
    """
    # Convert to Path objects
    protein_hits_path = Path(protein_hits_path)
    protein_metadata_path = Path(protein_metadata_path)
    chg_membership_path = Path(chg_membership_path)
    output_dir = Path(output_dir)

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read protein metadata
    print(f"Reading protein metadata from {protein_metadata_path}...")
    metadata = pd.read_csv(protein_metadata_path)

    # Create lookup dictionaries for faster access
    protein_to_genome = dict(zip(metadata["protein_id"], metadata["genome_accession"], strict=True))
    protein_to_taxon = dict(zip(metadata["protein_id"], metadata["taxon"], strict=True))

    # Read CHG membership
    print(f"Reading CHG membership from {chg_membership_path}...")
    chg_membership = pd.read_csv(chg_membership_path, sep="\t")

    # Create lookup dictionary for CHG IDs (we discard original_chg_id as requested)
    protein_to_chg = dict(zip(chg_membership["protein_id"], chg_membership["chg_id"], strict=True))

    # Read protein hits
    print(f"Reading protein hits from {protein_hits_path}...")
    hits = pd.read_csv(protein_hits_path, sep="\t")

    print(f"Total hits before filtering: {len(hits)}")

    # Add genome_accession for query and target
    hits["query_genome"] = hits["query"].map(protein_to_genome)
    hits["target_genome"] = hits["target"].map(protein_to_genome)

    # Add CHG IDs for query and target
    hits["chg"] = hits["query"].map(protein_to_chg)

    # Filter out hits where query and target have the same genome
    genome_hits = hits[hits["query_genome"] != hits["target_genome"]].copy()

    print(f"Hits after filtering same-genome hits: {len(genome_hits)}")

    # Add taxon information
    genome_hits["query_taxon"] = genome_hits["query"].map(protein_to_taxon)
    genome_hits["target_taxon"] = genome_hits["target"].map(protein_to_taxon)

    # Create normalized taxon pair (alphabetically sorted to avoid permutations)
    def create_taxon_pair(row):
        taxa = sorted([row["query_taxon"], row["target_taxon"]])
        return f"{taxa[0]}-{taxa[1]}"

    genome_hits["taxon_pair"] = genome_hits.apply(create_taxon_pair, axis=1)

    # Reorder columns for better readability
    columns_order = [
        "query",
        "target",
        "query_genome",
        "target_genome",
        "query_taxon",
        "target_taxon",
        "chg",
        "taxon_pair",
        "pident",
        "alnlen",
        "evalue",
        "bits",
    ]
    genome_hits = genome_hits[columns_order]

    # Save to output directory
    output_path = output_dir / "genome_hits.tsv"
    print(f"Saving genome hits to {output_path}...")
    genome_hits.to_csv(output_path, sep="\t", index=False)

    print(f"Done! Saved {len(genome_hits)} inter-genome hits.")
    print(f"Unique taxon pairs: {genome_hits['taxon_pair'].nunique()}")


def filter_max_pident_per_chg_taxon(genome_hits_path: str | Path, output_dir: str | Path) -> None:
    """
    Filter genome hits to keep only the maximum pident for each CHG-taxon pair combination.

    For each unique combination of CHG ID and taxon pair, this function keeps only
    the hit with the highest percent identity (pident).

    Parameters
    ----------
    genome_hits_path : str | Path
        Path to genome_hits.tsv containing columns:
        query, target, query_genome, target_genome, query_taxon, target_taxon,
        chg, taxon_pair, pident, alnlen, evalue, bits
    output_dir : str | Path
        Directory where filtered_genome_hits.tsv will be saved

    Returns
    -------
    None
        Saves filtered_genome_hits.tsv to the output directory
    """
    # Convert to Path objects
    genome_hits_path = Path(genome_hits_path)
    output_dir = Path(output_dir)

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read genome hits
    print(f"Reading genome hits from {genome_hits_path}...")
    genome_hits = pd.read_csv(genome_hits_path, sep="\t")

    print(f"Total hits before filtering: {len(genome_hits)}")
    print(
        f"Unique CHG-taxon pair combinations: {genome_hits.groupby(['chg', 'taxon_pair']).ngroups}"
    )

    # For each CHG-taxon pair combination, keep only the row with maximum pident
    # If there are ties (same max pident), keep the first one
    filtered_hits = genome_hits.loc[
        genome_hits.groupby(["chg", "taxon_pair"])["pident"].idxmax()
    ].copy()

    print(f"Hits after filtering (max pident per CHG-taxon pair): {len(filtered_hits)}")

    # Sort by CHG and taxon_pair for better readability
    filtered_hits = filtered_hits.sort_values(["chg", "taxon_pair"]).reset_index(drop=True)

    # Save to output directory
    output_path = output_dir / "filtered_genome_hits.tsv"
    print(f"Saving filtered genome hits to {output_path}...")
    filtered_hits.to_csv(output_path, sep="\t", index=False)

    print(f"Done! Saved {len(filtered_hits)} filtered hits.")
    print(f"Unique CHGs: {filtered_hits['chg'].nunique()}")
    print(f"Unique taxon pairs: {filtered_hits['taxon_pair'].nunique()}")

    # Show statistics about reduction
    reduction_pct = (1 - len(filtered_hits) / len(genome_hits)) * 100
    print(f"Reduction: {reduction_pct:.2f}% of original hits removed")


def bootstrap_confidence_intervals(  # noqa: PLR0913
    filtered_genome_hits_path: str | Path,
    output_path: str | Path,
    n_bootstrap: int = 10000,
    quantile: float = 0.95,
    confidence_level: float = 0.95,
    random_seed: int | None = 42,
) -> None:
    """
    Perform bootstrap analysis to determine confidence intervals for the 95th quantile
    of pident for each taxon pair and identify significantly high percent identities.

    This function:
    1. For each taxon pair, performs bootstrap resampling (default: 10,000 samples)
    2. Calculates the 95th percentile of pident for each bootstrap sample
    3. Determines the confidence interval for the 95th percentile
    4. Identifies hits with pident higher than the upper bound of the confidence interval

    Parameters
    ----------
    filtered_genome_hits_path : str | Path
        Path to filtered_genome_hits.tsv containing columns:
        query, target, query_genome, target_genome, query_taxon, target_taxon,
        chg, taxon_pair, pident, alnlen, evalue, bits
    output_path : str | Path
        Path where the results will be saved as a TSV file
    n_bootstrap : int, default=10000
        Number of bootstrap samples to generate
    quantile : float, default=0.95
        The quantile to calculate for each bootstrap sample (95th percentile)
    confidence_level : float, default=0.95
        Confidence level for the confidence interval (95%)
    random_seed : int | None, default=42
        Random seed for reproducibility. Set to None for non-reproducible results.

    Returns
    -------
    None
        Saves results to output_path with columns:
        - All original columns from filtered_genome_hits.tsv
        - q95_lower_bound: Lower bound of confidence interval for 95th percentile
        - q95_upper_bound: Upper bound of confidence interval for 95th percentile
        - is_significant: Boolean indicating if pident > upper bound
    """
    # Convert to Path objects
    filtered_genome_hits_path = Path(filtered_genome_hits_path)
    output_path = Path(output_path)

    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Set random seed for reproducibility
    if random_seed is not None:
        np.random.seed(random_seed)

    # Read filtered genome hits
    print(f"Reading filtered genome hits from {filtered_genome_hits_path}...")
    filtered_hits = pd.read_csv(filtered_genome_hits_path, sep="\t")

    print(f"Total hits: {len(filtered_hits)}")
    print(f"Unique taxon pairs: {filtered_hits['taxon_pair'].nunique()}")

    # Calculate alpha for confidence interval bounds
    alpha = 1 - confidence_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100

    print(f"\nPerforming bootstrap analysis with {n_bootstrap} samples...")
    print(f"Calculating {quantile*100}th percentile of pident for each taxon pair")
    print(
        f"Confidence interval: {confidence_level*100}% ({lower_percentile}% - {upper_percentile}%)"
    )

    # Store results for each taxon pair
    bootstrap_results = {}

    # Group by taxon pair and perform bootstrap
    for taxon_pair, group in filtered_hits.groupby("taxon_pair"):
        pident_values = group["pident"].values
        n_samples = len(pident_values)

        print(f"\n  Processing {taxon_pair}: {n_samples} hits")

        # Perform bootstrap sampling
        bootstrap_quantiles = np.zeros(n_bootstrap)
        for i in range(n_bootstrap):
            # Resample with replacement
            bootstrap_sample = np.random.choice(
                pident_values.astype(float), size=n_samples, replace=True
            )
            # Calculate the quantile for this bootstrap sample
            bootstrap_quantiles[i] = np.percentile(bootstrap_sample, quantile * 100)

        # Calculate confidence interval for the quantile
        ci_lower = np.percentile(bootstrap_quantiles, lower_percentile)
        ci_upper = np.percentile(bootstrap_quantiles, upper_percentile)

        bootstrap_results[taxon_pair] = {"q95_lower_bound": ci_lower, "q95_upper_bound": ci_upper}

        print(f"    95th percentile CI: [{ci_lower:.2f}, {ci_upper:.2f}]")

    # Add bootstrap results to the dataframe
    filtered_hits["q95_lower_bound"] = filtered_hits["taxon_pair"].map(
        lambda x: bootstrap_results[x]["q95_lower_bound"]
    )
    filtered_hits["q95_upper_bound"] = filtered_hits["taxon_pair"].map(
        lambda x: bootstrap_results[x]["q95_upper_bound"]
    )

    # Mark hits as significant if pident > upper bound of CI
    filtered_hits["is_significant"] = filtered_hits["pident"] > filtered_hits["q95_upper_bound"]

    # Filter to keep only significant hits
    significant_hits = filtered_hits[filtered_hits["is_significant"]].copy()

    # Save only significant results
    print(f"\nSaving significant hits to {output_path}...")
    significant_hits.to_csv(output_path, sep="\t", index=False)

    # Print summary statistics
    n_significant = filtered_hits["is_significant"].sum()
    pct_significant = (n_significant / len(filtered_hits)) * 100

    print("\n=== Summary ===")
    print(f"Total hits analyzed: {len(filtered_hits)}")
    print(f"Significant hits (pident > upper bound): {n_significant} ({pct_significant:.2f}%)")
    print(f"Significant hits saved to output file: {len(significant_hits)}")

    # Show breakdown by taxon pair
    print("\n=== Significant hits by taxon pair ===")
    sig_by_taxon = filtered_hits.groupby("taxon_pair").agg({"is_significant": ["sum", "count"]})
    sig_by_taxon.columns = ["significant", "total"]
    sig_by_taxon["percentage"] = (sig_by_taxon["significant"] / sig_by_taxon["total"]) * 100
    sig_by_taxon = sig_by_taxon.sort_values("significant", ascending=False)

    for taxon_pair, row in sig_by_taxon.iterrows():
        print(
            f"  {taxon_pair}: {int(row['significant'])}/{int(row['total'])} ({row['percentage']:.2f}%)"
        )

    print("\nDone!")


if __name__ == "__main__":
    # Example usage
    root = Path(__file__).parents[3] / "data" / "test"
    protein_hits_path = root / "03_alignment" / "hits_cleaned.tsv"
    protein_metadata_path = root / "02_preprocessing" / "protein_metadata.csv"
    chg_membership_path = root / "04_networks" / "reduced_chg_membership.tsv"
    output_dir = root / "07_lgt"

    # Step 1: Create genome hits
    create_genome_hits(
        protein_hits_path=protein_hits_path,
        protein_metadata_path=protein_metadata_path,
        chg_membership_path=chg_membership_path,
        output_dir=output_dir,
    )

    # Step 2: Filter to keep only max pident per CHG-taxon pair
    genome_hits_path = output_dir / "genome_hits.tsv"
    filter_max_pident_per_chg_taxon(genome_hits_path=genome_hits_path, output_dir=output_dir)

    # Step 3: Bootstrap confidence intervals for significant hits
    filtered_genome_hits_path = output_dir / "filtered_genome_hits.tsv"
    bootstrap_output_path = output_dir / "bootstrap_results.tsv"
    bootstrap_confidence_intervals(
        filtered_genome_hits_path=filtered_genome_hits_path,
        output_path=bootstrap_output_path,
        n_bootstrap=10000,
        quantile=0.95,
        confidence_level=0.95,
    )
