# Cancer-Selective Essential Gene Finder

Uses public CRISPR essentiality screens (BioGRID ORCS) to find genes that
look essential for cancer cell survival but are *not* flagged as essential
in normal-cell screens — i.e. candidate genes whose knockout might kill
cancer cells while sparing healthy tissue.

## What it does, step by step

1. **Filter the screen index** down to "pure" essentiality screens: gene
   knockout/inhibition screens, phenotype = proliferation/fitness, no added
   drug condition, scored with Bayes Factor (so scores are comparable across
   studies). This excludes drug-resistance, viral, and other unrelated
   screen types that happen to share the word "essential" in their
   rationale.
2. **Split screens into cancer vs. normal** cell lines using the `CELL_TYPE`
   column (matching Cancer/Leukemia/Carcinoma/Lymphoma/Melanoma/Glioma/Sarcoma).
3. **Read every screen file** and build a per-gene score map for each group
   (cancer, normal), keeping only rows where the Bayes Factor clears
   `BF_THRESHOLD` (default 5) — i.e. actual hits, not noise.
4. **Take the difference**: genes that were hits in cancer screens but never
   a hit in *any* normal-cell screen become the candidate list.
5. **Remove pan-essential genes** — genes required by almost every cell type
   (ribosomes, spliceosome, DNA replication machinery, etc.) aren't useful
   as selective targets even if they never showed up in our thin normal-cell
   sample. Filtered out using:
   - DepMap's `CRISPRInferredCommonEssentials` reference list
   - A secondary regex/name check for common housekeeping gene families
6. **Rank and export** the surviving candidates by two scores (see below),
   write them to CSV, and plot the top N of each.

## Two scoring methods

- **`average_score`** — mean Bayes Factor across the screens where the gene
  was a hit. Answers "how strongly essential is it, on average, when it
  does show up?"
- **`combined_score`** — `average_score * (screens hit / total screens)`.
  Answers "how strongly *and* how consistently essential is it across the
  whole cancer screen set?" Penalizes genes that only scored high in one or
  two screens.

## Results (this run)

See `outputs/cancer_specific_genes.csv` for the full ranked list, and
`outputs/plots/` for the top-10 bar charts under each scoring method.

**Read the results with caution — this is exploratory, not a validated
target list:**
- A gene passing every filter here means "well-supported by this dataset,"
  not "proven to kill cancer cells if knocked out." The Bayes Factor
  measures relative dropout in a pooled competition assay — a proxy for
  reduced fitness, which could mean cell death, growth arrest, or
  differentiation. Confirming actual cell death requires a follow-up assay.
- The "normal screen" side of the dataset is thin (~87 screens vs. ~275
  cancer screens), so "never a hit in normal screens" sometimes just means
  "wasn't tested enough," not "genuinely safe in normal cells."
- Any real candidate from this list should be cross-checked against the
  literature and DepMap's per-cell-line dependency data before being taken
  seriously as a target.

## Biological insights (short version)

The two output plots don't agree on a top gene, and that disagreement is
itself the interesting finding. `top_genes_avg_score.png` is topped by
**CRTC3**, a gene that only ever showed up as a hit in **1** screen — a
single extreme measurement, not a repeated one. `top_genes_combined_score.png`
(which rewards genes seen across *multiple* screens, not just a high single
score) instead puts **GART** on top, a gene backed by 8 independent screen
hits — a much more trustworthy signal.

Past GART and WRN, the combined-score chart flattens into a tight, low
cluster of genes that all look roughly equally "important" — that's a
noise floor, not a ranked signal. The root cause is dataset size: with
only 87 normal-cell screens, we can't reliably say a gene is "safe" in
normal cells just because it never showed up as a hit — it may simply
not have been tested enough. That's why this pipeline leans on DepMap's
much larger reference panel to catch pan-essential genes our own normal
screens were too sparse to flag.

**Full walkthrough of every exclusion step, with numbers, is in
[`BIOLOGICAL_INSIGHTS.md`](BIOLOGICAL_INSIGHTS.md).**

## How this could be used

- As a **first-pass filter** to shrink ~20,000 human genes down to a
  short, biologically-motivated candidate list worth manually researching.
- As a **teaching example** of the essentiality-screen analysis workflow:
  index filtering → per-screen aggregation → differential comparison →
  pan-essential exclusion.
- As a **starting point**, not an endpoint — real target nomination needs
  literature review, mutation-context awareness (e.g. is the target only
  a dependency in a specific cancer subtype), and experimental validation.

## Project structure

```
├── main.py                  # pipeline entry point — run this
├── config.py                # all paths and thresholds
├── src/
│   ├── data_loader.py       # screen index loading + filtering + cancer/normal split
│   ├── gene.py               # Gene class: accumulates scores across screens
│   ├── screen_processor.py   # reads screen files, builds gene_id -> Gene map
│   ├── filters.py            # cancer-specific / DepMap / housekeeping filters
│   └── visualize.py          # bar chart plots
├── outputs/
│   ├── cancer_specific_genes.csv
│   └── plots/
│       ├── top_genes_avg_score.png
│       └── top_genes_combined_score.png
├── screens_biogrid/          # data folder (not in git, see Setup below)
└── requirements.txt
```

## Setup

```bash
git clone <this-repo-url>
cd CRISPR
pip install -r requirements.txt
```

### Data (required, not included in the repo)

1. Download the full human BioGRID ORCS screen archive and extract every
   file into `screens_biogrid/`:
   https://downloads.thebiogrid.org/Download/BioGRID-ORCS/Latest-Release/BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens.tar.gz
   (~2.7GB, ~1900 files — this is why it isn't committed to git)
2. Download DepMap's common-essential gene list (`CRISPRInferredCommonEssentials.csv`)
   from the DepMap Data Portal (https://depmap.org/portal/download/) and
   place it in `screens_biogrid/` too. This step is optional — the pipeline
   still runs without it, just with a weaker pan-essential filter.

## Run

```bash
python main.py
```

Paths are resolved relative to the project root (via `config.py`), so this
works from any working directory, on any machine, as long as
`screens_biogrid/` is populated as described above.

## Configuration

All tunable values live in `config.py`:

| Parameter | Default | Meaning |
|---|---|---|
| `BF_THRESHOLD` | 5 | Minimum Bayes Factor to count a gene as a hit in one screen |
| `TOP_N` | 10 | How many top genes to include in the plots |

## Origin

This project started as exploratory analysis in `dryrun.ipynb` — the
notebook is kept in the repo as a record of the iterative process (including
dead ends, like discovering that raw score ranking mostly surfaces
housekeeping genes). `main.py` and `src/` are the cleaned-up, reproducible
version of that logic.

## Author

Shivani Patel — 2026
