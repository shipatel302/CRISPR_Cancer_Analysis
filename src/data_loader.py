"""Loads the BioGRID ORCS screen index and picks out the screens we want."""
import pandas as pd

# Cell type keywords used to split screens into cancer vs normal cell lines.
CANCER_KEYWORDS = "Cancer|Leukemia|Carcinoma|Lymphoma|Melanoma|Glioma|Sarcoma"


def load_screen_index(index_file):
    return pd.read_csv(index_file, sep="\t")


def filter_essentiality_screens(index_df):
    """Keep only screens that are:
    - about gene essentiality (not drug resistance, viral, etc.)
    - measuring cell proliferation
    - run without an added drug/condition (no confounding)
    - scored with Bayes Factor (so scores are comparable across screens)
    """
    essential = index_df[index_df["SCREEN_RATIONALE"].str.contains("essential", case=False, na=False)]
    essential = essential[essential["PHENOTYPE"].str.contains("proliferation", case=False, na=False)]
    essential = essential[essential["CONDITION_NAME"] == "-"]
    essential = essential[essential["SCORE.1_TYPE"] == "Bayes Factor"]
    return essential


def split_cancer_normal(essential_df):
    """Split filtered screens into cancer-cell-line screens and normal-cell-line screens."""
    cancer_ids = essential_df[essential_df["CELL_TYPE"].str.contains(CANCER_KEYWORDS, case=False, na=False)]["#SCREEN_ID"]
    normal_ids = essential_df[~essential_df["CELL_TYPE"].str.contains(CANCER_KEYWORDS, case=False, na=False)]["#SCREEN_ID"]
    return list(cancer_ids), list(normal_ids)


def screen_file_paths(data_dir, screen_ids):
    return [str(data_dir / f"BIOGRID-ORCS-SCREEN_{sid}-2.0.18.screen.tab.txt") for sid in screen_ids]
