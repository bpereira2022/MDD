# ============================================================
# QC Figures 
# ============================================================
library(DESeq2)
library(ggplot2)
library(vsn)
library(tidyverse)

# ── Load data ────────────────────────────────────────────────
metadata   <- read.csv("ColData_MDD.csv")
Count_Data <- read.csv("MDD_Final_Counts.csv",
                       row.names = 1, header = TRUE, check.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = Count_Data,
                              colData   = metadata,
                              design    = ~ 1)
dds <- dds[rowSums(counts(dds)) > 1, ]

# VST
vsd        <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_counts <- assay(vsd)

# Library sizes
raw_lib_sizes <- colSums(counts(dds))
vst_lib_sizes <- colSums(vst_counts)

# ── Shared theme constants ───────────────────────────────────
BASE_CEX_AXIS <- 1.1
BASE_CEX_LAB  <- 1.4
BASE_CEX_MAIN <- 1.6

# ── Cook's distance (needed for panels D and E) ──────────────
dds_fit <- estimateSizeFactors(dds)
dds_fit <- estimateDispersions(dds_fit)
dds_fit <- nbinomWaldTest(dds_fit)
cooks   <- assays(dds_fit)[["cooks"]]

# ── PCA data (needed for panel C) ───────────────────────────
pcaData    <- plotPCA(vsd, intgroup = colnames(metadata), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# ── Top 20 miRNAs (needed for panel F) ──────────────────────
expr_matrix     <- assay(vsd)
mean_expression <- rowMeans(expr_matrix)

top20 <- data.frame(
  miRNA     = rownames(expr_matrix),
  mean_expr = mean_expression
) %>%
  arrange(desc(mean_expr)) %>%
  slice_head(n = 20) %>%
  mutate(miRNA = factor(miRNA, levels = miRNA[order(mean_expr)]))

# ============================================================
# PANEL A — Raw library sizes
# ============================================================
png("QC_A_LibrarySizes_Raw.png", width = 3000, height = 2400, res = 300)
par(mar = c(9, 9, 4, 2),   # wide left margin for e+07 labels
    mgp = c(7, 0.8, 0))    # y-axis title pushed clear of tick labels
barplot(raw_lib_sizes,
        las       = 2,
        main      = "Library Sizes",
        ylab      = "Total Counts",
        col       = "grey60",
        border    = NA,
        cex.axis  = BASE_CEX_AXIS,
        cex.lab   = BASE_CEX_LAB,
        cex.main  = BASE_CEX_MAIN,
        cex.names = 0.75)
dev.off()

# ============================================================
# PANEL B — VST-transformed library sizes
# ============================================================
png("QC_B_LibrarySizes_VST.png", width = 3000, height = 2400, res = 300)
par(mar = c(9, 8, 4, 2),   # left margin for 5-digit tick labels
    mgp = c(6, 0.8, 0))    # y-axis title pushed clear of tick labels
barplot(vst_lib_sizes,
        las       = 2,
        main      = "VST-Transformed: Library Sizes",
        ylab      = "Sum of VST Values",
        col       = "grey60",
        border    = NA,
        cex.axis  = BASE_CEX_AXIS,
        cex.lab   = BASE_CEX_LAB,
        cex.main  = BASE_CEX_MAIN,
        cex.names = 0.75)
dev.off()

# ============================================================
# PANEL C — PCA plot
# ============================================================
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = metadata[[1]])) +
  geom_point(size = 3.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 11)
  ) +
  guides(color = guide_legend(ncol = 2))

ggsave("QC_C_PCA.png", plot = pca_plot, width = 10, height = 7, dpi = 300)

# ============================================================
# PANEL D — Cook's distance boxplot
# ============================================================
png("QC_D_CooksDistance.png", width = 3000, height = 1800, res = 300)
par(mar = c(9, 6, 4, 2),
    mgp = c(4, 0.8, 0))
boxplot(log10(cooks + 1),
        las      = 2,
        main     = "Cook's Distance (log10)",
        ylab     = "log10(Cook's Distance + 1)",
        col      = "grey80",
        cex.axis = BASE_CEX_AXIS,
        cex.lab  = BASE_CEX_LAB,
        cex.main = BASE_CEX_MAIN)
dev.off()

# ============================================================
# PANEL E — Dispersion estimates
# ============================================================
png("QC_E_Dispersion.png", width = 2400, height = 2000, res = 300)
par(mar     = c(5, 5, 4, 2),
    cex.axis = BASE_CEX_AXIS,
    cex.lab  = BASE_CEX_LAB,
    cex.main = BASE_CEX_MAIN)
plotDispEsts(dds_fit, main = "Dispersion Estimates")
dev.off()

# ============================================================
# PANEL F — Top 20 most expressed miRNAs
# ============================================================
top20_plot <- ggplot(top20, aes(miRNA, mean_expr)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Most Expressed miRNAs",
       x     = "miRNA",
       y     = "Mean Expression (VST Normalized)") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title  = element_text(size = 18, face = "bold"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 13)
  )

ggsave("QC_F_Top20miRNAs.png", plot = top20_plot,
       width = 8, height = 7, dpi = 300)

message("Done — 6 PNG files saved to your working directory.")
