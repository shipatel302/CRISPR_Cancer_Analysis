"""Reads individual screen files and builds a gene_id -> Gene map."""
import pandas as pd

from src.gene import Gene


def process_screens(file_list, threshold):
    """Read every screen file in file_list and record each gene's score,
    but only for rows whose Bayes Factor is >= threshold (i.e. actual hits).

    Returns (gene_map, error_screens, row_error_count, total_row_count).
    """
    gene_map = {}
    error_screens = []
    row_error_count = 0
    total_row_count = 0
    total_screens = len(file_list)

    for file in file_list:
        try:
            screen = pd.read_csv(file, sep="\t")
        except Exception:
            error_screens.append(file)
            continue

        for _, row in screen.iterrows():
            total_row_count += 1
            try:
                gene_id = int(row["IDENTIFIER_ID"])
                score = float(row["SCORE.1"])
            except (ValueError, TypeError):
                row_error_count += 1
                continue

            if score < threshold:
                continue

            if gene_id not in gene_map:
                gene_map[gene_id] = Gene(gene_id, symbol=row["OFFICIAL_SYMBOL"], total_screens=total_screens)

            gene_map[gene_id].add_score(score)
            gene_map[gene_id].add_repetition(row["#SCREEN_ID"])

    return gene_map, error_screens, row_error_count, total_row_count
