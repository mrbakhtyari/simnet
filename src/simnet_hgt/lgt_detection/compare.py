"""
Module for comparing and merging HGT analysis results.
"""

from pathlib import Path

import pandas as pd


def merge_bootstrap_with_hgt_scores(
    bootstrap_results_path: str | Path, hgt_scores_path: str | Path, output_path: str | Path
) -> None:
    """
    Merge bootstrap results with HGT scores and report minimum HGT score.

    This function:
    1. Reads bootstrap results (significant hits) and HGT scores
    2. Matches hits based on protein pairs (query/target with protein_A/protein_B)
    3. Adds HGT_score column to bootstrap results
    4. Reports the minimum HGT_score found
    5. Saves the merged results

    Parameters
    ----------
    bootstrap_results_path : str | Path
        Path to bootstrap_results.tsv containing columns:
        query, target, chg, and other bootstrap analysis columns
    hgt_scores_path : str | Path
        Path to hgt_scores.tsv containing columns:
        protein_A, protein_B, CHG, HGT_score, and other HGT analysis columns
    output_path : str | Path
        Path where the merged results will be saved as a TSV file

    Returns
    -------
    None
        Saves merged results to output_path with HGT_score added
    """
    # Convert to Path objects
    bootstrap_results_path = Path(bootstrap_results_path)
    hgt_scores_path = Path(hgt_scores_path)
    output_path = Path(output_path)

    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Read input files
    print(f"Reading bootstrap results from {bootstrap_results_path}...")
    bootstrap = pd.read_csv(bootstrap_results_path, sep="\t")

    print(f"Reading HGT scores from {hgt_scores_path}...")
    hgt_scores = pd.read_csv(hgt_scores_path, sep="\t")

    print(f"\nBootstrap results: {len(bootstrap)} rows")
    print(f"HGT scores: {len(hgt_scores)} rows")

    # Create a lookup dictionary for HGT scores
    # Key: (protein_A, protein_B) - matching query/target from bootstrap
    print("\nCreating HGT score lookup...")
    hgt_lookup = {}

    for _, row in hgt_scores.iterrows():
        # Store both directions since query/target order might differ
        key1 = (row["protein_A"], row["protein_B"])
        key2 = (row["protein_B"], row["protein_A"])
        hgt_lookup[key1] = row["HGT_score"]
        hgt_lookup[key2] = row["HGT_score"]

    # Match bootstrap results with HGT scores
    print("Matching bootstrap results with HGT scores...")
    hgt_score_values = []
    matches_found = 0

    for _, row in bootstrap.iterrows():
        key = (row["query"], row["target"])
        hgt_score = hgt_lookup.get(key)

        if hgt_score is not None:
            matches_found += 1

        hgt_score_values.append(hgt_score)

    # Add HGT_score column to bootstrap results
    bootstrap["HGT_score"] = hgt_score_values

    # Report statistics
    print("\nMatching statistics:")
    print(
        f"  Matches found: {matches_found}/{len(bootstrap)} ({matches_found/len(bootstrap)*100:.2f}%)"
    )
    print(f"  Unmatched: {len(bootstrap) - matches_found}")

    # Filter out rows without HGT scores and report
    bootstrap_with_scores = bootstrap[bootstrap["HGT_score"].notna()].copy()

    if len(bootstrap_with_scores) > 0:
        min_hgt_score = bootstrap_with_scores["HGT_score"].min()
        max_hgt_score = bootstrap_with_scores["HGT_score"].max()
        mean_hgt_score = bootstrap_with_scores["HGT_score"].mean()

        print("\nHGT Score statistics:")
        print(f"  Minimum HGT_score: {min_hgt_score:.6f}")
        print(f"  Maximum HGT_score: {max_hgt_score:.6f}")
        print(f"  Mean HGT_score: {mean_hgt_score:.6f}")

        # Show some examples
        print("\n=== Examples of merged results ===")
        print("Top 5 by HGT_score:")
        top_5 = bootstrap_with_scores.nlargest(5, "HGT_score")[
            ["chg", "taxon_pair", "pident", "HGT_score", "is_significant"]
        ]
        print(top_5.to_string(index=False))

        print("\nBottom 5 by HGT_score:")
        bottom_5 = bootstrap_with_scores.nsmallest(5, "HGT_score")[
            ["chg", "taxon_pair", "pident", "HGT_score", "is_significant"]
        ]
        print(bottom_5.to_string(index=False))
    else:
        print("\nWARNING: No matches found between bootstrap results and HGT scores!")

    # Save results
    print(f"\nSaving merged results to {output_path}...")
    bootstrap_with_scores.to_csv(output_path, sep="\t", index=False)

    print(f"\nDone! Saved {len(bootstrap_with_scores)} rows with HGT scores.")


if __name__ == "__main__":
    # Example usage
    root = Path(__file__).parents[3] / "data" / "test"
    bootstrap_results_path = root / "07_lgt" / "bootstrap_results.tsv"
    hgt_scores_path = root / "06_hgt" / "hgt_scores.tsv"
    output_path = root / "07_lgt" / "significant_hgt_candidates.tsv"

    merge_bootstrap_with_hgt_scores(
        bootstrap_results_path=bootstrap_results_path,
        hgt_scores_path=hgt_scores_path,
        output_path=output_path,
    )
