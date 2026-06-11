# MDD EV miRNA Analysis Repository

This repository contains the processing and analysis outputs for extracellular vesicle (EV) miRNA sequencing data from a Major Depressive Disorder (MDD) blood plasma study across exercise-related timepoints.

---

## Repository Contents

- `/MDD Processing`
  Linux shell pipeline for:
  1. Adapter trimming (`bbduk.sh`)
  2. Alignment to human miRNA reference (`bowtie`)
  3. Count generation from aligned reads (`samtools`)

- `/FastQC`
  Raw and trimmed read quality-control reports (`fastqc`, `multiqc`).

- `/DEG Analysis Results`
  Differential expression outputs by timepoint:
  - `DE_mid.csv`
  - `DE_post.csv`
  - `DE_75min.csv`
  - `Table1_DE_miRNAs.csv`
  - `expressed_miRNAs.txt`
  - `Shared_DEGs_table.csv`

- `/Pathway Analysis Results`
  Pathway enrichment summaries and plots:
  - `Pathway Results.pdf`
  - `bubble_top_all.csv`
  - `category_counts.csv`
  - `panther_categorized.csv`

- `/Target Results`
  Target-gene and overlap files:
  - `target_mid.tsv`
  - `target_post.tsv`
  - `target_75.tsv`
  - `panther_mid.txt`
  - `panther_post.txt`
  - `panther_75.txt`
  - `target_sets.csv`
  - `target_overlap_summary.csv`

---

## Data Files in Repository Root

- `MDD_Final_Counts.csv` — miRNA count matrix (rows = miRNAs, columns = samples)
- `ColData_MDD.csv` — sample metadata (sample IDs, treatment timepoints, subject IDs)

---

## R Scripts

### `figure_colors.R`
Defines all shared color palettes used across figures. Source this file at the top of any analysis script.

- **Timepoint colors** (Okabe-Ito color-blind-safe): PRE = black, MID = bluish-green (`#009E73`), POST = blue (`#0072B2`), 75 min POST = vermillion (`#D55E00`)
- **Pretty timepoint labels** for legends and facets (e.g., `"MID"` → `"Mid-AAE"`)
- **DEG pattern annotation colors** for Figure 4A row annotation (single and multi-timepoint patterns)
- **`divergent_palette(n)`** — symmetric blue-white-red color ramp for log2FC heatmaps (Figures 4A, 4B)
- **`divergent_breaks(cap, n)`** — generates symmetric breaks for consistent heatmap scaling
- **`jaccard_palette(n)`** — white-to-blue sequential palette for Jaccard similarity heatmap (Figure 5B)

---

### `Quality_Check.R`
Exploratory QC workflow. Run this before differential expression analysis to inspect data quality and flag potential outliers.

**Inputs:** `MDD_Final_Counts.csv`, `ColData_MDD.csv`

**Steps:**
1. Build a DESeqDataSet (design `~ 1`, no model) and filter low-count miRNAs (`rowSums > 1`)
2. **Library size check** — raw count barplot and VST-transformed barplot, plus a side-by-side comparison
3. **Variance-stabilizing transformation (VST)** — blind VST for downstream QC plots
4. **PCA** — plots PC1 vs PC2, colored by first metadata column
5. **Sample distance heatmap** — hierarchical clustering of Euclidean distances in VST space (`pheatmap`, Blues palette)
6. **Dispersion estimates** — `plotDispEsts()` to verify negative binomial fit
7. **Mean-SD plot** — `meanSdPlot()` from `vsn` to check heteroscedasticity
8. **Cook's distance** — boxplot (log10) to flag high-influence samples
9. **Top 20 expressed miRNAs** — bar chart of mean VST expression
10. **Isolation Forest outlier detection** — `isotree` package; flags samples in top 5% anomaly score
11. **Outlier conclusion** — sample `352_E1` flagged by isolation forest but not confirmed by Cook's distance or PCA; all samples retained

**Key packages:** `DESeq2`, `pheatmap`, `RColorBrewer`, `ggplot2`, `vsn`, `isotree`, `ggrepel`

---

### `QC_figures.R`
Produces publication-ready PNG panels from the QC analysis. Outputs 6 files to the working directory.

**Inputs:** `MDD_Final_Counts.csv`, `ColData_MDD.csv`

**Output files:**

| File | Description |
|------|-------------|
| `QC_A_LibrarySizes_Raw.png` | Raw count library sizes (barplot) |
| `QC_B_LibrarySizes_VST.png` | VST-transformed library sizes (barplot) |
| `QC_C_PCA.png` | PCA plot, colored by first metadata column |
| `QC_D_CooksDistance.png` | Cook's distance boxplot (log10) |
| `QC_E_Dispersion.png` | Dispersion estimates plot |
| `QC_F_Top20miRNAs.png` | Top 20 most expressed miRNAs (mean VST) |

**Key packages:** `DESeq2`, `ggplot2`, `vsn`, `tidyverse`

---

### `Master_Analysis.R`
Main differential expression and figure generation pipeline. Run after QC is complete.

**Inputs:** `MDD_Final_Counts.csv`, `ColData_MDD.csv`, `figure_colors.R`, `target_sets.csv` (optional, for Figure 5)

**All outputs are saved to a `/results` subdirectory.**

**Key parameters (set at top of script):**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_COUNT` | 10 | Minimum normalized counts for low-count filter |
| `MIN_SAMPLES` | 3 | Minimum number of samples meeting `MIN_COUNT` |
| `FDR_SIG` | 0.05 | Significance threshold |
| `FDR_TREND` | 0.10 | Trending threshold |

**Steps and outputs:**

1. **Load & align data** — verifies that all count columns have matching ColData rows; reorders ColData to match count column order
2. **DESeq2 model** — design `~ ID + Treatment`; PRE as reference; low-count filter applied; non-converging genes (~2%) dropped
3. **Three contrasts** — MID vs PRE, POST vs PRE, 75 min POST vs PRE
4. **Table 1** (`Table1_DE_miRNAs.csv`) — all significant and trending DE miRNAs across all three contrasts, ordered by FDR; also writes `DE_mid.csv`, `DE_post.csv`, `DE_75min.csv`
5. **Figure 3A** (`Figure_3A_PCA.tiff`) — faceted PCA with 95% ellipses; three panels (PRE vs each post-exercise timepoint)
6. **Figure 3B** (`Figure_3B_Volcano_mid/post/75.tiff`) — EnhancedVolcano plots for each contrast; top 12 significant miRNAs labeled
7. **Per-subject paired log2FC matrix** — computed from normalized counts, used by Figures 4A and 4B
8. **Figure 4A** (`Figure_4A_Heatmap_all.tiff`) — heatmap of all significantly DE miRNAs; rows annotated by DEG pattern (which timepoints), columns by timepoint
9. **Figure 4B** (`Figure_4B_Heatmap_shared.tiff`) — heatmap of miRNAs trending (FDR < 0.10) in ≥2 timepoints; rows ordered by mean log2FC
10. **Figure 4C** (`Figure_4C_Venn_DE.tiff`) — Venn diagram of trending DE miRNAs across the three contrasts; `Shared_DEGs_table.csv` written with per-timepoint log2FC and FDR
11. **Figure 5** (`Figure_5A_Venn_targets.tiff`, `Figure_5B_Jaccard.tiff`) — target gene overlap Venn and Jaccard similarity heatmap; skipped if `target_sets.csv` is not present
12. **Figure 7 (partial)** — normalized count boxplot for `hsa-miR-1252-3p` across all timepoints; saved as PDF and TIFF

**Key packages:** `DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `ggrepel`, `pheatmap`, `RColorBrewer`, `EnhancedVolcano`, `VennDiagram`, `viridis`

---

### `Figures_6_7.R`
Generates pathway enrichment figures and the focused miR-513a-3p panels. Run after `Master_Analysis.R` and after PANTHER pathway results are available.

**Inputs:** `figure_colors.R`, `bubble_top10_all.csv` (optional), `pathway_categories.csv` (optional), `ColData_Final.csv`, `MDD_Final_Counts.csv`

**Steps and outputs:**

| Figure | File | Description |
|--------|------|-------------|
| 6A | `Figure_6A_bubble.tiff` | Bubble plot: top 10 enriched GO terms per timepoint; size = fold enrichment, color = −log10(Bonferroni); skipped if `bubble_top10_all.csv` not found |
| 6B | `Figure_6B_pies.tiff` | Faceted pie charts of pathway categories per timepoint; skipped if `pathway_categories.csv` not found |
| 7A | `Figure_7A_table.tiff` / `Figure_7A_gene_table.csv` | Table of shared dendrite-morphogenesis target genes with directional regulation across timepoints; rendered as TIFF if `ggpubr` is installed |
| 7B | `Figure_7B_513a_log2FC.tiff` | Bar chart of DESeq2 log2FC ± lfcSE for `hsa-miR-513a-3p` vs Pre-AAE at each timepoint |
| 7C | `Figure_7C_513a_bubble.tiff` | Bubble plot: top 10 enriched GO pathways for `hsa-miR-513a-3p` targets |

**Key packages:** `ggplot2`, `dplyr`, `tidyr`, `tibble`, `readr`, `viridis`, `DESeq2`, `ggpubr` (optional)

---

## Recommended Run Order

```
1. figure_colors.R       # sourced automatically by other scripts; no need to run standalone
2. Quality_Check.R       # exploratory QC
3. QC_figures.R          # publication QC panels
4. Master_Analysis.R     # DE analysis and Figures 3–5, partial Figure 7
5. Figures_6_7.R         # pathway and miR-513a figures (requires PANTHER outputs)
```

---

## Notes

- The processing script in `/MDD Processing` contains local path placeholders and should be updated for your environment before running.
- `Figure_6A` and `Figure_6B` require PANTHER output files (`bubble_top10_all.csv`, `pathway_categories.csv`) to be present in the working directory; those are in the folder /Pathway Analysis Results.
- `Figure_5` requires `target_sets.csv` (columns: `mid`, `post`, `min75`) exported from DIANA-microT (score ≥ 0.8); in /Target Results.
- `Figure_7B` in `Figures_6_7.R` requires `ColData_Final.csv`.
