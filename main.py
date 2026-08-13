"""
Find genes that look essential for cancer cell survival but are not
flagged as essential in normal-cell CRISPR screens.

Data source: BioGRID ORCS human CRISPR screens (see README for download link).

Run with:
    python main.py
"""
import pandas as pd

# importing config like index directory path and screens path
import config


# IMport utility functions
#
from src.data_loader import (
    load_screen_index,
    filter_essentiality_screens,
    split_cancer_normal,
    screen_file_paths,
)

from src.screen_processor import process_screens

# filters to remove essential genes for normal cells
# cancer_specific_genes removes normal essential genes from cancer essential genes from screens
# load_depmap_common_essentials gives set of pan essential genes from CRISPR
# exclude_ids function to exclude normal gene from cancerous gene
# exclude_housekeeping removes normal known pan essential genes from cancer genes
from src.filters import (
    cancer_specific_genes,
    load_depmap_common_essentials,
    exclude_ids,
    exclude_housekeeping,
)
from src.visualize import plot_top_genes


def main():
    #check if output directories are present else crete them
    config.OUTPUT_DIR.mkdir(exist_ok=True)
    config.PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: narrow the screen index down to pure essentiality screens,
    # then split into cancer-cell-line screens vs normal-cell-line screens.
    index_df = load_screen_index(config.INDEX_FILE)
    essential = filter_essentiality_screens(index_df)

    # AFter filtering essential screens we want to divide them in Cancer and Normal Screen IDs
    cancer_ids, normal_ids = split_cancer_normal(essential)
    print(f"Cancer essentiality screens: {len(cancer_ids)}")
    print(f"Normal essentiality screens: {len(normal_ids)}")

    # Create list of files based upon screen id f"/BIOGRID-ORCS-SCREEN_{sid}-2.0.18.screen.tab.txt
    cancer_files = screen_file_paths(config.DATA_DIR, cancer_ids)
    normal_files = screen_file_paths(config.DATA_DIR, normal_ids)

    # Creating map of Cancerous Gene and Normal gene reading diffrent screen files 

    gene_map_cancer, _, cancer_row_errors, cancer_rows = process_screens(cancer_files, config.BF_THRESHOLD)
    gene_map_normal, _, normal_row_errors, normal_rows = process_screens(normal_files, config.BF_THRESHOLD)
    print(f"Cancer screens: {cancer_rows} rows read, {cancer_row_errors} skipped (bad gene id / score)")
    print(f"Normal screens: {normal_rows} rows read, {normal_row_errors} skipped (bad gene id / score)")

    # Step 3: keep genes essential in cancer screens but never a hit in normal-cell screens.
    candidates = cancer_specific_genes(gene_map_cancer, gene_map_normal)
    print(f"Cancer-specific candidates (before pan-essential filtering): {len(candidates)}")

    # Step 4: remove genes known to be pan-essential  CRISPR  (needed by almost every cell type,
    # cancer or not) using DepMap's reference list plus a housekeeping-gene name check.
    if config.DEPMAP_ESSENTIALS_FILE.exists():
        depmap_ids = load_depmap_common_essentials(config.DEPMAP_ESSENTIALS_FILE)
        candidates = exclude_ids(candidates, depmap_ids)
    else:
        print(f"Warning: {config.DEPMAP_ESSENTIALS_FILE} not found, skipping DepMap filter")

    candidates = exclude_housekeeping(candidates)
    print(f"Cancer-selective candidates (final): {len(candidates)}")

    # Step 5: write the ranked candidate list to CSV.
    rows = [
        {
            "symbol": g.symbol,
            "entrez_id": g.identifier_id,
            "average_score": g.average_score(),
            "combined_score": g.combined_score(),
            "num_screens_hit": g.num_repetitions(),
        }
        for g in candidates
    ]
    results_df = pd.DataFrame(rows).sort_values("average_score", ascending=False)
    results_csv = config.OUTPUT_DIR / "cancer_specific_genes.csv"
    results_df.to_csv(results_csv, index=False)
    print(f"Saved {len(results_df)} genes to {results_csv}")

    # Step 6: plot the top genes by each scoring method.
    plot_top_genes(
        candidates,
        lambda g: g.average_score(),
        "Top candidate genes by average Bayes Factor",
        "Average Bayes Factor",
        config.PLOTS_DIR / "top_genes_avg_score.png",
        top_n=config.TOP_N,
    )
    plot_top_genes(
        candidates,
        lambda g: g.combined_score(),
        "Top candidate genes by combined score (avg score x consistency)",
        "Combined score",
        config.PLOTS_DIR / "top_genes_combined_score.png",
        top_n=config.TOP_N,
    )
    print(f"Saved plots to {config.PLOTS_DIR}")


if __name__ == "__main__":
    main()
