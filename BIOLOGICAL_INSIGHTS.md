# Biological Insights

## Objective

The goal of this project is to identify **differentially essential genes**
between cancer and normal cells using public CRISPR knockout screen data —
in plain terms, genes a cancer cell cannot survive without, but a normal
cell can. That's the working definition of a good drug target: knocking
it out should collapse the tumor while leaving healthy tissue largely
unharmed.

Scientifically, this means comparing gene essentiality (measured by Bayes
Factor scores from BAGEL analysis of CRISPR dropout screens) across two
populations of screens — cancer cell lines vs. normal cell lines  and
isolating genes with a strong essentiality signal in the cancer group but
no corresponding signal in the normal group. In practice, that comparison
is only as good as the "normal" side of the data, which is where most of
the interesting problems in this project came from.


## The problem:

I was looking for genes that cancer cells need to survive but normal cells don't so if you knocked that gene out, the cancer would die and healthy cells would be fine.

The hard part is proving a gene is actually safe to remove from normal cells. To know that for sure, I need to test it on enough normal cells. If I don't, I can't tell "this gene really isn't needed" apart from "I just didn't check enough." And if I get that wrong, my whole list isn't really cancer-specific  it's just genes every cell needs to live.

## What I excluded, step by step

**1,952 screens → 362 essentiality screens.**
Most of the screens in BioGRID ORCS have nothing to do with essentiality —
plenty are drug resistance screens, viral screens, differentiation screens.
I kept only the ones tagged "essential" in rationale, measuring
proliferation, run without a drug, and scored with Bayes Factor, so I'm
comparing apples to apples: every remaining screen is asking the same
question — *does losing this gene stop the cell from growing?*

**362 screens → 275 cancer / 87 normal.**
I split these by cell type. And this is where I ran into my biggest
limitation: 87 normal screens just isn't a lot. There are around 20,000
human genes — 87 screens can't possibly give me confident coverage of all
of them. If a gene never shows up as a hit in those 87 screens, that could
mean it's genuinely safe in normal cells, or it could just mean it was
never tested enough. I couldn't tell those two apart from my own data
alone.

**Scoring:** I counted a gene as a "hit" in a screen if its Bayes Factor
was ≥ 5. Across everything, I read about 4.93M rows from cancer screens
(992 skipped for bad data) and 1.56M rows from normal screens (895
skipped).

**Cancer hits minus normal hits → 1,063 candidates.**
Any gene that showed up as a hit in a cancer screen and never showed up as
a hit in any of the 87 normal screens made it into this list. But given
how thin the normal side is, I knew a good chunk of these 1,063 genes
were probably just genes I hadn't tested enough on the normal side — not
genes that are actually cancer-specific.

**So I brought in an outside reference: DepMap's common-essential gene
list.** DepMap has run CRISPR knockouts across 1,000+ cell lines — way
more coverage than my 87 screens — and published which genes come out
essential almost everywhere (ribosomes, splicing machinery, DNA
replication, that kind of thing). I used that list to catch the
housekeeping genes my own normal-screen sample was too small to flag on
its own. This one step removed most of the obvious false positives.

**Then a smaller cleanup pass on gene names.** I added a regex to catch a
few more housekeeping gene families by naming pattern (`RPL*`, `RPS*`,
`*ARS`, `PSM*`, histones, etc.) that slipped past the DepMap list. I'll
be upfront that this is a blunt tool — it only catches genes whose name
happens to reveal their function, so it's a patch, not a real filter.

**What was left: 985 candidate genes.**

## What the plots are actually showing

I made two rankings, and they don't agree — which turned out to be the
most useful thing I found.

The **average score chart** ranks genes by their mean Bayes Factor. The
top gene, CRTC3, scored 156 — but when I checked, it was only ever a hit
in **1** screen. Same story for a few others near the top (TRIM8, TTC9C,
BCAS4, HABP2). A single extreme number from one screen could be a real
finding, or it could just be noise — with only one data point, I have no
way to tell the difference. So this chart is really ranking "biggest
single measurement," not "most reliable gene."

The **combined score chart** multiplies the average score by how
consistently the gene shows up across screens. That changes the winner:
**GART** takes the top spot here, even though its average score (98) is
lower than CRTC3's (156) — because GART was a hit in 8 separate screens,
not 1. Eight screens agreeing is a much stronger signal than one screen
with a big number.

Here's the part I think matters most: look at how the combined-score chart
drops off. GART (2.85) and WRN (1.42) stand somewhat apart, and then
everything after that — FH, TIMM17A, ALDH18A1, CTNNB1, and the rest — sits
in a tight, low band (roughly 0.2 to 0.7) that barely separates from each
other. That flat cluster looks like noise to me, not a real ranking. If
my filtering had actually isolated cancer-specific genes cleanly, I'd
expect to see a handful of genes standing clearly above a baseline.
Instead I get two candidates that look real, and then a long tail I can't
distinguish from chance.

## What I'd actually say about this

- **GART and WRN** are the two candidates I'd defend from this run — both
  have multiple independent screens backing them up, not just one lucky
  measurement. WRN also happens to have outside support: it's already a
  real synthetic-lethal drug target being pursued for MSI-high cancers,
  which lines up with what the data is pointing at.
- Everything else near the top of either list is resting on 1–2 screen
  hits, which isn't enough for me to trust on its own.
- The core issue isn't my filtering logic — it's that I only had 87
  normal-cell screens to work with. That's not enough to say "safe in
  normal cells" with real confidence, which is exactly why I had to lean
  on DepMap's much bigger reference panel to catch what my own data
  couldn't.
- If I wanted to trust this pipeline without leaning on an outside
  reference, I'd need a lot more normal-cell screen data — enough that
  "never a hit in normal screens" actually means something, instead of
  "wasn't tested enough times to know."
