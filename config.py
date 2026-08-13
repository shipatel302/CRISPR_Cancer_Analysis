"""All paths and thresholds for the pipeline, in one place."""
from pathlib import Path

# Paths are relative to this file, not the current working directory,
# so the project runs the same no matter where it's launched from.
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "screens_biogrid"
OUTPUT_DIR = BASE_DIR / "outputs"
PLOTS_DIR = OUTPUT_DIR / "plots"

INDEX_FILE = DATA_DIR / "BIOGRID-ORCS-SCREEN_INDEX-2.0.18.index.tab.txt"
DEPMAP_ESSENTIALS_FILE = DATA_DIR / "CRISPRInferredCommonEssentials.csv"

# Minimum Bayes Factor for a gene to count as "essential" in a single screen.
BF_THRESHOLD = 5

# How many top genes to write to the plots.
TOP_N = 10
