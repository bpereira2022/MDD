# MDD EV miRNA Analysis Repository

This repository contains the processing and analysis outputs for extracellular vesicle (EV) miRNA sequencing data from a Major Depressive Disorder (MDD) blood plasma study across exercise-related timepoints.

## Repository contents

- `/MDD Processing`  
  Linux shell pipeline for:
  1. Adapter trimming (`bbduk.sh`)
  2. Alignment to human miRNA reference (`bowtie`)
  3. Count generation from aligned reads (`samtools`)

- `/FastQC`  
  Raw and trimmed read quality-control reports (`fastqc`, `multiqc`).

- `/Quality_Check.R`  
  R workflow for:
  - QC diagnostics (library size checks, PCA, sample-distance heatmap, dispersion/mean-SD plots)
  - Outlier detection (Cook’s distance + isolation forest)
  - Top expressed miRNA visualization

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

## Data files in repository root

- `MDD_Final_Counts.csv` — miRNA count matrix
- `ColData_MDD.csv` — sample metadata

## Notes

- The processing script in `/MDD Processing` contains local path placeholders and should be updated for your environment before running.
- The R workflow currently references `ColData_Final_Anhodenia.csv`; update this filename/path to your available metadata file if needed.
