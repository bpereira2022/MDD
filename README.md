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
  FastQC/MultiQC workflow script for generating raw and trimmed read quality-control reports.

- `/DEG Analysis Results`
  Differential expression outputs by timepoint:
  - `DE_mid.csv`
  - `DE_post.csv`
  - `DE_75min.csv`
  - `Table1_DE_miRNAs.csv`
  - `expressed_miRNAs.txt`

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

### `hemolysis_platelet_QC.R`
Post-hoc QC workflow for hemolysis and platelet contamination signals using normalized counts from the DESeq2 object.

**Inputs:** `ddsClean` object from the main DE pipeline (`Master_Analysis.R`)

**Output files (saved to `/results`):**

| File | Description |
|------|-------------|
| `Supplementary_Figure2_Hemolysis_QC.tiff` / `.pdf` | Hemolysis index boxplot + PCA colored by hemolysis index |
| `Supplementary_Table_HemolysisIndex.csv` | Per-sample hemolysis index values |
| `Supplementary_Table_PlateletMiRNA_FriedmanTests.csv` | Per-miRNA Friedman test results for platelet-enriched miRNAs |
| `Supplementary_Figure3_Platelet_QC.tiff` / `.pdf` | Faceted platelet-enriched miRNA expression boxplots across timepoints |

**Key checks:**
1. Computes hemolysis index from canonical numerator/denominator miRNAs and tests timepoint effects with a Friedman test
2. Flags hemolysis-threshold context in visualization (`log2(7)` reference line)
3. Evaluates platelet-enriched miRNAs across timepoints with per-miRNA Friedman tests (BH-adjusted p-values)
4. Highlights platelet miRNAs that are also DE in the main study results

**Key packages:** `DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `patchwork`

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
10. **Figure 4C** (`Figure_4C_Venn_DE.tiff`) — Venn diagram of trending DE miRNAs across the three contrasts; script also exports a per-timepoint overlap summary table
11. **Figure 5** (`Figure_5A_Venn_targets.tiff`, `Figure_5B_Jaccard.tiff`) — target gene overlap Venn and Jaccard similarity heatmap; skipped if `target_sets.csv` is not present
12. **Figure 7 (partial)** — normalized count boxplot for `hsa-miR-1252-3p` across all timepoints; saved as PDF and TIFF

**Key packages:** `DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `ggrepel`, `pheatmap`, `RColorBrewer`, `EnhancedVolcano`, `VennDiagram`, `viridis`

---

## Recommended Run Order

```
1. figure_colors.R       # sourced automatically by other scripts; no need to run standalone
2. Quality_Check.R       # exploratory QC
3. QC_figures.R          # publication QC panels
4. Master_Analysis.R     # DE analysis and Figures 3–5, partial Figure 7
5. hemolysis_platelet_QC.R  # post-hoc hemolysis and platelet contamination QC (after ddsClean is available)
```

---

## Notes

- The processing script in `/MDD Processing` contains local path placeholders and should be updated for your environment before running.
- `Figure_5` in `Master_Analysis.R` requires `target_sets.csv` (columns: `mid`, `post`, `min75`) exported from DIANA-microT (score ≥ 0.8); the script will skip Figure 5 with an informative message if the file is missing.
