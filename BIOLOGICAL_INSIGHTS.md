# Biological Insights

This is the "why" behind the pipeline — what we're actually filtering out,
why we had to borrow an outside gene list to do it properly, and how to
read the plots without over-trusting them.

## The problem in plain words

We want genes that a cancer cell can't live without, but a normal cell can.
That second half — "but a normal cell can" — is the hard part, because it
means we need to know what normal cells *do* need, and rule those genes out.

## Step by step: what got excluded and why

1. **Screen index (1,952 screens total) → 362 essentiality screens.**
   Most CRISPR screens in BioGRID ORCS aren't about essentiality at all —
   they're drug resistance screens, viral infection screens, differentiation
   screens, etc. We kept only screens rationale-tagged "essential",
   phenotype "proliferation", no drug condition, and scored with Bayes
   Factor, so every remaining screen is asking the same question in a
   comparable way: *does losing this gene stop the cell growing?*

2. **362 essentiality screens → 275 cancer / 87 normal.**
   Split by matching cancer-related keywords in `CELL_TYPE`. This is the
   first place things get shaky: **87 normal screens is a small
   sample.** There are roughly 20,000 human genes. 87 screens, even at
   full coverage each, is nowhere near enough independent looks to say
   with confidence "this gene is never needed by a normal cell." Absence
   of evidence in 87 screens is weak evidence of absence.

3. **Per-screen gene scoring.** Every gene-row with Bayes Factor ≥ 5 in a
   screen counts as a "hit" in that screen. Across all screens: ~4.93M
   cancer-screen rows read (992 skipped — bad IDs/scores), ~1.56M
   normal-screen rows read (895 skipped).

4. **Cancer hits minus normal hits → 1,063 candidates.**
   Any gene that was a hit in at least one cancer screen and *never* a hit
   in any of the 87 normal screens survives this step. Given point 2, a
   lot of these 1,063 "survivors" are really just genes that never got
   properly tested on the normal side, not genes that are truly safe in
   normal cells.

5. **This is exactly why we brought in DepMap's `CRISPRInferredCommonEssentials`
   list.** DepMap ran CRISPR knockouts across 1,000+ cell lines — a much
   bigger and more systematic sample than our 87 normal screens — and
   published the genes that come out essential almost everywhere
   (ribosomes, spliceosome, DNA replication machinery, etc.). We use that
   list as an external stand-in for "things our own normal-screen sample
   was too thin to catch." This step alone removes most of the
   false-positive housekeeping genes that our own comparison missed.

6. **A second, smaller symbol-pattern filter catches a few more.**
   Gene family naming patterns (`RPL*`, `RPS*`, `*ARS`, `PSM*`, spliceosome
   components, histones, etc.) catch a handful of additional pan-essential
   genes that weren't in the DepMap list for one reason or another. This
   step is a blunt instrument — it can only catch genes whose *name*
   reveals their function, so it's a cleanup pass, not a real filter.

7. **Final list: 985 candidates**, ranked by two scores (see README).

## Reading the plots

**`top_genes_avg_score.png`** — ranks genes purely by their mean Bayes
Factor. Look at the `num_screens_hit` column in the CSV next to this
chart's top genes: **CRTC3 (top bar, score 156) was only ever a hit in
1 screen.** So was TRIM8, TTC9C, BCAS4, HABP2, INO80C. A single high
score from one screen can be a real biological signal, or it can be one
noisy measurement — with n=1 you can't tell the difference. This chart
is mostly ranking by "how extreme was the single best measurement," not
"how reliable is this gene."

**`top_genes_combined_score.png`** — multiplies average score by how
often the gene showed up across all screens (consistency). This is why
the ranking changes: **GART jumps to #1** here even though its average
score (98) was lower than CRTC3's (156) — GART was a hit in 8 separate
screens, not 1. That's a much stronger signal: 8 independent screens
agreeing beats one screen with an extreme number.

Look at the shape of the combined-score chart: there's a real drop after
GART (2.85) and WRN (1.42), then everything from FH onward sits in a
tight, low cluster (roughly 0.2–0.7) that barely separates from each
other. That flat cluster is the tell — **it looks like a noise floor, not
a ranked list of real signals.** If the pipeline were cleanly separating
"genuinely cancer-selective" from "not," you'd expect a handful of genes
clearly standing above a baseline. Instead we get two-ish standouts
(GART, WRN) and then a long tail that's statistically indistinguishable
from chance.

## What this actually tells us

- **GART** and **WRN** are the only two candidates in this run with
  multi-screen support behind them, which makes them the most defensible
  results from this pipeline — not proof, but the best-supported entries.
  (WRN in particular has real external validation: it's an actively
  pursued synthetic-lethal target in MSI-high cancers — see the notebook
  discussion.)
- Everything else in the top 10–20 is built on 1–2 screen hits, which is
  too thin a base to trust on its own.
- The deeper issue isn't the filtering logic — it's the **input data
  size**. With only 87 normal screens, we structurally cannot build a
  reliable "essential in normal cells" reference from our own data alone;
  that's why an external reference (DepMap) had to fill the gap, and even
  that only helps with genes DepMap happened to flag as common-essential.
- **Bottom line: to trust this list without leaning on outside references,
  we'd need a much larger and more diverse normal-cell screen collection**
  — enough independent screens per gene that "never a hit in normal
  screens" actually means something, instead of "wasn't tested enough
  times to know."
