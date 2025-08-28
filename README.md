MDD Processing: Code used in Linux system to 1) remove adaptor sequence (bbduk.sh), 2) alingment to the human miRNA reference genome (bowtie1), and generation of count data (samtools).

FastQC: Quality check of raw and trimed reads (fastqc v0.12.1, multiqc v1.14). 

Final_Code_July15: R Markdown workflow for quality control, outlier detection, differential expression analysis, and visualization of extracellular vesicle (EV) miRNA sequencing data from a study on MDD and blood plasma. The pipeline is designed to analyze changes in miRNA expression across multiple post-exercise timepoints:
1) Quality Control: Performs sample clustering, principal component analysis (PCA), mean-variance diagnostics, and outlier identification using Cook’s distance and isolation forest.
2) Differential Expression Analysis:
   - Creates a DESeq2 dataset with experimental design (subject ID and treatment/timepoint).
   - Normalizes counts, removes genes failing model convergence, and performs differential expression testing between timepoints (PRE, MID, POST, 75 min POST).
   - Outputs significant results to CSV files.
3) Visualization:
   - Generates PCA plots to visualize sample separation by treatment/timepoint.
   - Creates volcano plots highlighting significant miRNAs.
   - Plots heatmaps for selected and all significant miRNAs, showing log2 fold change across subjects and timepoints.
   - Bubble plots visualize pathway enrichment results for top GO terms and specific miRNA targets.
4) Jaccard Similarity and Gene Overlap:
   - Calculates Jaccard similarity indices for miRNA target gene sets between timepoints.
   - Summarizes and visualizes the overlap of predicted target genes for selected biological processes.

Data Requirements:
- Count matrix (`MDD_Final_Counts.csv`) (provided)
- Sample metadata (`ColData_Final_Anhodenia.csv`) (provided)
- Pathway enrichment results and gene sets for similarity analysis (See Supplementary File with pathway enrichment results).
