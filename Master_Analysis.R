## Packages
# install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "EnhancedVolcano"))
# install.packages(c("ggplot2","dplyr","tidyr","tibble","ggrepel","readr",
#                    "pheatmap","RColorBrewer","VennDiagram","viridis","ggpubr"))

## Open packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(VennDiagram)
  library(viridis)
  library(grid)
})

## ===========================================================================
## Setting Up (:
## ===========================================================================
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
count     <- dplyr::count
slice     <- dplyr::slice

# Shared color-blind-safe palette 
source("figure_colors.R")

# Output 
dir.create("results", showWarnings = FALSE)
out <- function(f) file.path("results", f)

#   low-count filter = keep miRNAs with >= MIN_COUNT normalized counts in >= MIN_SAMPLES samples.
MIN_COUNT   <- 10
MIN_SAMPLES <- 3

#   cut offs
FDR_SIG     <- 0.05   # significant
FDR_TREND   <- 0.10   # trending

## ===========================================================================
## Load Data
## ===========================================================================
metadata   <- read.csv("ColData_MDD.csv", fileEncoding = "UTF-8-BOM")
Count_Data <- read.csv("MDD_Final_Counts.csv", row.names = 1,
                       header = TRUE, check.names = FALSE)

cat("Loaded:", nrow(Count_Data), "miRNAs x", ncol(Count_Data), "samples\n")


## ===========================================================================
## Align Samples
## ===========================================================================

stopifnot("Sample_ID" %in% colnames(metadata))
rownames(metadata) <- metadata$Sample_ID

# Are all count columns represented in ColData?
missing <- setdiff(colnames(Count_Data), rownames(metadata))
if (length(missing) > 0)
  stop("These count columns have no ColData row: ", paste(missing, collapse = ", "))

# Reorder ColData rows to match count columns exactly
metadata <- metadata[colnames(Count_Data), , drop = FALSE]

# Refuse to continue unless perfectly aligned
if (!all(rownames(metadata) == colnames(Count_Data)))
  stop("ALIGNMENT FAILED - colData rows do not match count columns. Stopping.")

frac_already <- mean(metadata$Sample_ID == colnames(Count_Data))  # post-reorder = 1
cat("Sample alignment check passed. (Reordered ColData to count columns.)\n")


## ===========================================================================
## Build DESeq2 object, filter, fit
## ===========================================================================
metadata$ID        <- factor(metadata$ID)
metadata$Treatment <- factor(metadata$Treatment,
                             levels = c("PRE", "MID", "POST", "75 min POST"))

dds <- DESeqDataSetFromMatrix(countData = Count_Data,
                              colData   = metadata,
                              design    = ~ ID + Treatment)
dds <- estimateSizeFactors(dds)

# Low-count filter (normalized counts), then PRE as reference
keep <- rowSums(counts(dds, normalized = TRUE) >= MIN_COUNT) >= MIN_SAMPLES
dds  <- dds[keep, ]
dds$Treatment <- relevel(dds$Treatment, ref = "PRE")
dds  <- DESeq(dds)

# convergence issues check
cat("miRNAs passing filter (>=", MIN_COUNT, "in >=", MIN_SAMPLES, "samples):",
    nrow(dds), "\n")
cat("Convergence:\n"); print(table(mcols(dds)$betaConv))

# Drop genes that failed to converge (~2%)
ddsClean <- dds[which(mcols(dds)$betaConv), ]

cat("Samples per timepoint:\n"); print(table(colData(ddsClean)$Treatment))

# Export the full set of expressed miRNAs 
writeLines(rownames(ddsClean), out("expressed_miRNAs.txt"))
cat("Wrote", out("expressed_miRNAs.txt"), "(", nrow(ddsClean), "expressed miRNAs )\n")


## ===========================================================================
## Contrasts (each post-exercise timepoint vs PRE)
## ===========================================================================
res_mid  <- results(ddsClean, contrast = c("Treatment", "MID",         "PRE"))
res_post <- results(ddsClean, contrast = c("Treatment", "POST",        "PRE"))
res_75   <- results(ddsClean, contrast = c("Treatment", "75 min POST", "PRE"))

# Report DE counts 
de_count <- function(res, a) sum(na.omit(res)$padj < a)
cat(sprintf("\nSignificant (FDR < %.2f):  mid=%d  post=%d  75min=%d\n",
            FDR_SIG, de_count(res_mid, FDR_SIG), de_count(res_post, FDR_SIG),
            de_count(res_75, FDR_SIG)))
cat(sprintf("Trending   (FDR < %.2f):  mid=%d  post=%d  75min=%d\n",
            FDR_TREND, de_count(res_mid, FDR_TREND), de_count(res_post, FDR_TREND),
            de_count(res_75, FDR_TREND)))


## ===========================================================================
## TABLE 1  (DE + trending miRNAs, all three contrasts)
## ===========================================================================
## One tidy CSV with a Timepoint column, ordered by FDR within each block, including everything up to the trending threshold (FDR < 0.10).
table1_block <- function(res, label) {
  as.data.frame(res) %>%
    rownames_to_column("miRNA") %>%
    filter(!is.na(padj), padj < FDR_TREND) %>%
    transmute(Timepoint = label,
              miRNA,
              baseMean       = round(baseMean, 1),
              log2FoldChange = round(log2FoldChange, 2),
              padj           = signif(padj, 2)) %>%
    arrange(padj)
}
table1 <- bind_rows(
  table1_block(res_mid,  "Mid-AAE vs Pre-AAE"),
  table1_block(res_post, "Post-AAE vs Pre-AAE"),
  table1_block(res_75,   "75 min Post-AAE vs Pre-AAE"))
write_csv(table1, out("Table1_DE_miRNAs.csv"))

# Also write each contrast separately (handy for downstream DIANA input)
write_csv(table1_block(res_mid,  "Mid-AAE vs Pre-AAE"),  out("DE_mid.csv"))
write_csv(table1_block(res_post, "Post-AAE vs Pre-AAE"), out("DE_post.csv"))
write_csv(table1_block(res_75,   "75 min Post-AAE vs Pre-AAE"), out("DE_75min.csv"))
cat("Wrote Table 1 ->", out("Table1_DE_miRNAs.csv"), "\n")


## ===========================================================================
## FIGURE 3A  -  PCA (Pre vs each timepoint, variance ellipses)
## ===========================================================================
vsd        <- varianceStabilizingTransformation(ddsClean, blind = FALSE)
pca_data   <- plotPCA(vsd, intgroup = "Treatment", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

mk <- function(keep, lab) pca_data %>% filter(Treatment %in% keep) %>% mutate(Facet = lab)
subset_data <- bind_rows(
  mk(c("PRE", "MID"),         "Pre-AAE and Mid-AAE"),
  mk(c("PRE", "POST"),        "Pre-AAE and Post-AAE"),
  mk(c("PRE", "75 min POST"), "Pre-AAE and 75 minutes post-AAE"))
subset_data$Facet <- factor(subset_data$Facet,
  levels = c("Pre-AAE and Mid-AAE", "Pre-AAE and Post-AAE",
             "Pre-AAE and 75 minutes post-AAE"))

fig3a <- ggplot(subset_data, aes(PC1, PC2, color = Treatment, fill = Treatment)) +
  stat_ellipse(type = "norm", level = 0.95, geom = "polygon",
               alpha = 0.18, linewidth = 0.6) +     # 0.68 for ~1-SD bubbles
  geom_point(size = 2.5) +
  scale_color_manual(values = tp_colors, labels = tp_labels) +
  scale_fill_manual(values = tp_colors, guide = "none") +
  facet_wrap(~ Facet, ncol = 1, scales = "free") +
  labs(title = "PCA of Serum EV miRNA Expression: Pre and AAE Time Points",
       x = paste0("PC1: ", percentVar[1], "% Variance"),
       y = paste0("PC2: ", percentVar[2], "% Variance"),
       color = "Timepoint") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 11))

ggsave(out("Figure_3A_PCA.tiff"), fig3a, width = 5, height = 9,
       dpi = 300, compression = "lzw")


## ===========================================================================
## FIGURE 3B  -  Volcano plots 
## ===========================================================================
volcano_one <- function(res, title, subtitle,
                        xrange = c(-7, 7), yrange = c(0, 4.3), n_label = 12) {
  res_df <- as.data.frame(na.omit(res))
  sig    <- rownames(res_df)[res_df$padj < FDR_SIG]
  if (length(sig) > n_label) {
    ord <- order(res_df[sig, "padj"])
    sig <- sig[ord][seq_len(n_label)]
  }
  EnhancedVolcano(
    res_df, lab = rownames(res_df), x = "log2FoldChange", y = "padj",
    selectLab = sig, title = title, subtitle = subtitle, caption = NULL,
    pCutoff = FDR_SIG, FCcutoff = 0.1, pointSize = 1.0, labSize = 4.5,
    colAlpha = 0.6,
    legendLabels = c("NS", "Log2FC", "FDR", "FDR + Log2FC"),
    legendPosition = "top",
    drawConnectors = TRUE, widthConnectors = 0.4,
    min.segment.length = 0, max.overlaps = Inf,
    col = c("grey70", "#56B4E9", "#009E73", "#D55E00"),
    ylab = bquote(~-Log[10] ~ italic(FDR)),
    xlab = bquote(~Log[2] ~ "Fold Change")) +
    scale_x_continuous(breaks = seq(xrange[1], xrange[2], 2)) +
    coord_cartesian(xlim = xrange, ylim = yrange, clip = "off") +
    theme(plot.margin = margin(8, 30, 8, 8))
}
v_mid  <- volcano_one(res_mid,  "Mid-AAE vs. Pre-AAE",
                      "How does EV miRNA expression change during exercise?")
v_post <- volcano_one(res_post, "Post-AAE vs. Pre-AAE",
                      "What is the immediate response once exercise ends?")
v_75   <- volcano_one(res_75,   "75 Minutes Post-AAE vs. Pre-AAE",
                      "Do exercise-induced changes persist into recovery?")
ggsave(out("Figure_3B_Volcano_mid.tiff"),  v_mid,  width = 7.5, height = 6.5, dpi = 500, compression = "lzw")
ggsave(out("Figure_3B_Volcano_post.tiff"), v_post, width = 7.5, height = 6.5, dpi = 500, compression = "lzw")
ggsave(out("Figure_3B_Volcano_75.tiff"),   v_75,   width = 7.5, height = 6.5, dpi = 500, compression = "lzw")


## ===========================================================================
## Per-subject paired log2FC (shared by Fig 4A, 4B)
## ===========================================================================
norm_counts <- counts(ddsClean, normalized = TRUE) %>% as.data.frame()
norm_counts$miRNA <- rownames(norm_counts)
meta_df <- as.data.frame(colData(ddsClean)) %>%
  rownames_to_column("sample") %>% select(sample, ID, Treatment)

long_avg <- norm_counts %>%
  pivot_longer(-miRNA, names_to = "sample", values_to = "norm_count") %>%
  left_join(meta_df, by = "sample") %>%
  group_by(miRNA, ID, Treatment) %>%
  summarise(mean_count = mean(norm_count), .groups = "drop")

wide_counts <- long_avg %>% pivot_wider(names_from = Treatment, values_from = mean_count)

get_log2fc <- function(df, tp) df %>%
  filter(!is.na(PRE), !is.na(.data[[tp]])) %>%
  mutate(log2FC = log2((.data[[tp]] + 1) / (PRE + 1))) %>%
  transmute(miRNA, ID, log2FC, Timepoint = tp)

log2fc_all <- bind_rows(get_log2fc(wide_counts, "MID"),
                        get_log2fc(wide_counts, "POST"),
                        get_log2fc(wide_counts, "75 min POST"))


## ===========================================================================
## FIGURE 4A  -  Heatmap of all significant DE miRNAs
## ===========================================================================
## Rows ordered by mean log2FC: most upregulated at top, most down at bottom.
get_sig <- function(res, label) as.data.frame(res) %>%
  rownames_to_column("miRNA") %>% filter(!is.na(padj), padj < FDR_SIG) %>%
  mutate(DE_timepoint = label)
sig_all <- bind_rows(get_sig(res_mid, "MID"), get_sig(res_post, "POST"),
                     get_sig(res_75, "75MIN"))

de_matrix <- sig_all %>% select(miRNA, DE_timepoint) %>% distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = DE_timepoint, values_from = present, values_fill = 0)
for (cc in c("MID","POST","75MIN")) if (!cc %in% names(de_matrix)) de_matrix[[cc]] <- 0
de_matrix <- de_matrix %>%
  mutate(DE_pattern = case_when(
    MID==1 & POST==0 & `75MIN`==0 ~ "Mid-AAE Only",
    MID==0 & POST==1 & `75MIN`==0 ~ "Post-AAE Only",
    MID==1 & POST==1 & `75MIN`==0 ~ "Mid & Post-AAE",
    MID==0 & POST==0 & `75MIN`==1 ~ "75min-AAE Only",
    MID==1 & POST==0 & `75MIN`==1 ~ "Mid & 75min-AAE",
    MID==0 & POST==1 & `75MIN`==1 ~ "Post & 75min-AAE",
    MID==1 & POST==1 & `75MIN`==1 ~ "All timepoints",
    TRUE ~ "Other")) %>%
  arrange(factor(DE_pattern, levels = c(
    "Mid-AAE Only","Post-AAE Only","Mid & Post-AAE","75min-AAE Only",
    "Post & 75min-AAE","Mid & 75min-AAE","All timepoints")))

lfc4a <- log2fc_all %>% filter(miRNA %in% de_matrix$miRNA) %>%
  mutate(sample_id = paste(ID, Timepoint, sep = "_")) %>%
  select(miRNA, sample_id, log2FC) %>%
  pivot_wider(names_from = sample_id, values_from = log2FC) %>%
  column_to_rownames("miRNA") %>% as.matrix()
lfc4a <- lfc4a[de_matrix$miRNA[de_matrix$miRNA %in% rownames(lfc4a)], , drop = FALSE]

ca <- data.frame(sample_id = colnames(lfc4a)) %>%
  separate(sample_id, into = c("ID","Timepoint"), sep = "_(?=[^_]+$)", remove = FALSE)
ord <- order(match(ca$Timepoint, c("MID","POST","75 min POST")))
lfc4a <- lfc4a[, ord]; ca <- ca[ord, ]
annotation_col <- data.frame(
  Timepoint = factor(ca$Timepoint, levels = c("MID","POST","75 min POST")),
  row.names = ca$sample_id)
annotation_row <- data.frame(
  `DEG Pattern` = factor(de_matrix$DE_pattern),
  row.names = de_matrix$miRNA, check.names = FALSE)
ann_colors <- list(
  Timepoint = c("MID" = unname(tp_colors["MID"]),
                "POST" = unname(tp_colors["POST"]),
                "75 min POST" = unname(tp_colors["75 min POST"])),
  `DEG Pattern` = deg_pattern_colors)

tiff(out("Figure_4A_Heatmap_all.tiff"), width = 11, height = 9, units = "in",
     res = 300, compression = "lzw")
pheatmap(lfc4a, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = ann_colors, color = divergent_palette(100),
         breaks = divergent_breaks(cap = 10, n = 100),
         show_colnames = FALSE, fontsize_row = 7, na_col = "grey90",
         main = "Significantly DE Serum EV miRNAs Over Time in AAE")
dev.off()


## ===========================================================================
## FIGURE 4B  -  Heatmap of shared DE miRNAs (>=2 timepoints, trending)
## ===========================================================================
## Shared = trending (FDR < 0.10) in at least two contrasts.
trend <- function(res) rownames(na.omit(res))[na.omit(res)$padj < FDR_TREND]
tally <- table(c(trend(res_mid), trend(res_post), trend(res_75)))
shared <- names(tally)[tally >= 2]
cat("Shared DE miRNAs (trending, >=2 timepoints):", length(shared), "\n")

if (length(shared) >= 1) {
  m4b <- log2fc_all %>% filter(miRNA %in% shared) %>%
    mutate(Timepoint = recode(Timepoint, "75 min POST" = "75minPOST"),
           sample_id = paste(ID, Timepoint, sep = "_")) %>%
    select(miRNA, sample_id, log2FC) %>%
    pivot_wider(names_from = sample_id, values_from = log2FC) %>%
    column_to_rownames("miRNA") %>% as.matrix()
  cb <- data.frame(sample_id = colnames(m4b)) %>%
    separate(sample_id, into = c("ID","Timepoint"), sep = "_(?=[^_]+$)", remove = FALSE)
  o <- order(match(cb$Timepoint, c("MID","POST","75minPOST")))
  m4b <- m4b[, o]; cb <- cb[o, ]
  m4b <- m4b[names(sort(rowMeans(m4b, na.rm = TRUE), decreasing = TRUE)), , drop = FALSE]
  acol <- data.frame(Timepoint = factor(cb$Timepoint, levels = c("MID","POST","75minPOST")),
                     row.names = cb$sample_id)
  acolors <- list(Timepoint = c("MID" = unname(tp_colors["MID"]),
                                "POST" = unname(tp_colors["POST"]),
                                "75minPOST" = unname(tp_colors["75 min POST"])))
  tiff(out("Figure_4B_Heatmap_shared.tiff"), width = 11, height = 4.5, units = "in",
       res = 300, compression = "lzw")
  pheatmap(m4b, cluster_rows = FALSE, cluster_cols = FALSE,
           color = divergent_palette(100), breaks = divergent_breaks(cap = 8, n = 100),
           na_col = "grey90", annotation_col = acol, annotation_colors = acolors,
           show_colnames = FALSE, fontsize_row = 10,
           main = "Shared Differentially Expressed Serum EV miRNAs Over Time in AAE")
  dev.off()
}


## ===========================================================================
## FIGURE 4C  -  Venn of shared DE miRNAs (trending, FDR < 0.10)
## ===========================================================================
venn_de <- list("Mid-AAE" = trend(res_mid),
                "Post-AAE" = trend(res_post),
                "75 Minutes Post-AAE" = trend(res_75))
venn.diagram(venn_de, filename = out("Figure_4C_Venn_DE.tiff"),
             imagetype = "tiff", height = 2200, width = 2200, resolution = 300,
             compression = "lzw",
             fill = c(unname(tp_colors["MID"]), unname(tp_colors["POST"]),
                      unname(tp_colors["75 min POST"])),
             alpha = 0.35, cex = 1.4, cat.cex = 1.2, lwd = 2,
             main = "Shared Differentially Expressed miRNAs")
# venn.diagram writes a log file; remove it quietly
unlink(list.files("results", pattern = "VennDiagram.*\\.log$", full.names = TRUE))

## --- Shared-DEG names table ----------- ##
## For every miRNA trending (FDR<0.10) in >=2 timepoints, record its log2FC and FDR at each timepoint and which timepoints it is shared across.
trend_tab <- function(res, label) as.data.frame(res) %>%
  rownames_to_column("miRNA") %>%
  filter(!is.na(padj), padj < FDR_TREND) %>%
  transmute(miRNA, !!paste0(label, "_log2FC") := round(log2FoldChange, 2),
            !!paste0(label, "_FDR") := signif(padj, 2))
shared_tbl <- Reduce(function(a, b) full_join(a, b, by = "miRNA"),
                     list(trend_tab(res_mid,  "Mid"),
                          trend_tab(res_post, "Post"),
                          trend_tab(res_75,   "x75min"))) %>%
  rowwise() %>%
  mutate(n_timepoints = sum(!is.na(c_across(ends_with("_FDR"))))) %>%
  ungroup() %>%
  filter(n_timepoints >= 2) %>%
  arrange(desc(n_timepoints), miRNA)
# label which timepoints each miRNA is shared across (row-wise)
shared_tbl$shared_in <- apply(
  !is.na(shared_tbl[, c("Mid_FDR","Post_FDR","x75min_FDR")]), 1,
  function(v) paste(c("Mid","Post","75min")[v], collapse = ", "))
write_csv(shared_tbl, out("Shared_DEGs_table.csv"))
cat("Wrote", out("Shared_DEGs_table.csv"), "(", nrow(shared_tbl),
    "miRNAs shared across >=2 timepoints )\n")


## ===========================================================================
## FIGURE 5  -  Target overlap (needs DIANA target_sets.csv); in github
## ===========================================================================
if (file.exists("target_sets.csv")) {
  ts <- read.csv("target_sets.csv", check.names = FALSE)
  gene_sets <- list("Mid-AAE"  = unique(na.omit(ts$mid)),
                    "Post-AAE" = unique(na.omit(ts$post)),
                    "75 Minutes Post-AAE" = unique(na.omit(ts$min75)))

  venn.diagram(gene_sets, filename = out("Figure_5A_Venn_targets.tiff"),
               imagetype = "tiff", height = 2200, width = 2200, resolution = 300,
               compression = "lzw",
               fill = c(unname(tp_colors["MID"]), unname(tp_colors["POST"]),
                        unname(tp_colors["75 min POST"])),
               alpha = 0.40, cex = 1.5, cat.cex = 1.2, lwd = 2,
               main = "Shared Predicted Targets Across Time")
  unlink(list.files("results", pattern = "VennDiagram.*\\.log$", full.names = TRUE))

  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  nm <- names(gene_sets); jac <- matrix(NA, 3, 3, dimnames = list(nm, nm))
  for (i in nm) for (j in nm) jac[i, j] <- jaccard(gene_sets[[i]], gene_sets[[j]])
  tiff(out("Figure_5B_Jaccard.tiff"), width = 6, height = 5, units = "in",
       res = 300, compression = "lzw")
  pheatmap(jac, cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.2f", number_color = "black",
           fontsize_number = 12, color = jaccard_palette(100),
           main = "Jaccard Similarity of Targets Across Time in AAE")
  dev.off()
} else {
  message("[Fig 5] target_sets.csv not found - skipping. Export DIANA-microT ",
          "targets (score >= 0.8) as columns mid, post, min75.")
}


## ===========================================================================
## FIGURE 7  -  hsa-miR-1252-3p (7A)
## ===========================================================================

FOCUS_MIR <- "hsa-miR-1252-3p"

## 7A: focus-miRNA effect = DESeq2 model estimate +/- lfcSE across timepoints
library(ggplot2)
library(dplyr)

# Your color palette
tp_colors <- c(
  "PRE" = "black",
  "MID" = "#009E73",
  "POST" = "#0072B2",
  "75 min POST" = "#D55E00"
)

if (FOCUS_MIR %in% rownames(ddsClean)) {
  
  # Extract normalized counts for this miRNA across all samples
  norm_counts <- counts(ddsClean, normalized = TRUE)[FOCUS_MIR, ]
  
  # Create dataframe with sample info and counts
  count_data <- data.frame(
    counts = norm_counts,
    sample = names(norm_counts)
  ) %>%
    left_join(
      data.frame(sample = rownames(colData(ddsClean)),
                 Treatment = ddsClean$Treatment,
                 ID = ddsClean$ID),
      by = "sample"
    )
  
  # Map treatment to timepoint labels
  count_data$Timepoint <- recode(count_data$Treatment,
                                 "PRE" = "PRE",
                                 "MID" = "MID",
                                 "POST" = "POST",
                                 "75 min POST" = "75 min POST")
  
  count_data$Timepoint <- factor(count_data$Timepoint, 
                                 levels = c("PRE", "MID", "POST", "75 min POST"))
  
  # Create figure with boxplot + points (log scale with pseudocount)
  fig_counts <- ggplot(count_data, aes(x = Timepoint, y = counts + 0.5, fill = Timepoint)) +
    geom_boxplot(alpha = 0.7, color = "black", linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 2, color = "grey20") +
    scale_y_log10(expand = c(0.1, 0)) +
    scale_fill_manual(values = tp_colors, guide = "none") +
    labs(title = paste0(FOCUS_MIR, " Normalized Counts Across Timepoints"),
         x = "Timepoint",
         y = "Normalized Counts + 0.5 (log10)") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.3)
    )
  
  ggsave(paste0(FOCUS_MIR, "_normalized_counts.pdf"), fig_counts, 
         width = 6, height = 5, dpi = 300)
  ggsave(paste0(FOCUS_MIR, "_normalized_counts.tiff"), fig_counts, 
         width = 6, height = 5, dpi = 300, compression = "lzw")
  
  print(fig_counts)
  cat("\n✓ Figure saved!\n")
  
} else {
  message("[Counts Plot] ", FOCUS_MIR, " not in filtered set - skipping.")
}
