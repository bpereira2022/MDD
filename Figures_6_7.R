suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(viridis)
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

source("figure_colors.R")   # tp_colors, viridis_option, divergent_*, etc.

theme_pub <- theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 9),
        plot.title  = element_text(size = 11, face = "bold"),
        legend.position = "right")

## ===========================================================================
## Figure 6
## ===========================================================================

## ---- 6A. Bubble plot: top-10 pathways across timepoints ------------------ ##
if (file.exists("bubble_top10_all.csv")) {

  bubble_data <- read_csv("bubble_top10_all.csv", show_col_types = FALSE)

  # Order GO terms by timepoint, then significance (most significant at top)
  bubble_data$`GO Biological Process Complete Term` <- factor(
    bubble_data$`GO Biological Process Complete Term`,
    levels = rev(unique(bubble_data %>%
      arrange(`Time Point`, Bonferroni) %>%
      pull(`GO Biological Process Complete Term`))))

  bubble_data$`Time Point` <- factor(bubble_data$`Time Point`,
    levels = c("Mid-AAE", "Post-AAE", "75 minutes Post-AAE"))

  fig6a <- ggplot(bubble_data,
                  aes(x = `Time Point`, y = `GO Biological Process Complete Term`)) +
    geom_point(aes(size = `Fold Enrichment`, color = -log10(Bonferroni))) +
    scale_color_viridis_c(option = viridis_option, name = expression(-log[10]~"(Bonferroni)")) +
    scale_size_continuous(name = "Fold Enrichment", range = c(2, 8)) +
    labs(title = "Top Enriched Pathways of Target Genes at Each Timepoint",
         x = "Time Point", y = "Enriched GO Term") +
    theme_pub

  ggsave("Figure_6A_bubble.tiff", fig6a, width = 8.5, height = 8,
         dpi = 300, compression = "lzw")
} else {
  message("[6A] bubble_top10_all.csv not found - skipping. ",
          "Export your PANTHER top-10 table with columns: 'Time Point', ",
          "'GO Biological Process Complete Term', 'Fold Enrichment', 'Bonferroni', 'Count'.")
}

## ---- 6B. Pie charts: pathway categories per timepoint -------------------- ##
category_colors <- c(
  "Neurodevelopment and Maintenance" = "#009E73",
  "Cardiac Development"              = "#D55E00",
  "RNA Processing"                  = "#0072B2",
  "Signaling and Regulation"        = "#E69F00",
  "Muscle and Growth"               = "#CC79A7",
  "Migration and Morphogenesis"     = "#56B4E9",
  "Secretion and Import"            = "#F0E442"
)

if (file.exists("pathway_categories.csv")) {

  cats <- read_csv("pathway_categories.csv", show_col_types = FALSE)
  cats$`Time Point` <- factor(cats$`Time Point`,
    levels = c("Mid-AAE", "Post-AAE", "75 minutes Post-AAE"))
  cats$Category <- factor(cats$Category, levels = names(category_colors))

  fig6b <- ggplot(cats, aes(x = "", y = Count, fill = Category)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar(theta = "y") +
    facet_wrap(~ `Time Point`) +
    scale_fill_manual(values = category_colors, name = "Pathway Category") +
    labs(title = "Pathway Categories by Timepoint") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.position = "right",
          strip.text = element_text(face = "bold"))

  ggsave("Figure_6B_pies.tiff", fig6b, width = 10, height = 4,
         dpi = 300, compression = "lzw")
} else {
  message("[6B] pathway_categories.csv not found - skipping. ",
          "Provide columns: 'Time Point', 'Category', 'Count' ",
          "(Category one of: ", paste(names(category_colors), collapse = ", "), ").")
}

## ===========================================================================
## Figure 7
## ===========================================================================

## ---- 7A. Table of shared dendrite-morphogenesis target genes ------------- ##
fig7a_tbl <- tibble::tribble(
  ~Gene,     ~Function,                                            ~`Mid-AAE`, ~`Post-AAE`, ~`75 min Post-AAE`,
  "NRP1",    "Angiogenesis, axon guidance, neuronal development",   "Up",       "Down",      "Down",
  "PTEN",    "PI3K/Akt signaling; synaptic plasticity & survival",  "Up",       "Up",        "Up",
  "MEF2A",   "Synaptic plasticity; distinguishes MDD",              "Up",       "Up",        "Up",
  "SLITRK5", "Synaptic function; dysregulated in MDD",              "Up",       "Up",        "Up",
  "TMEM106B","Lysosomal pathways; dysregulated in MDD",             "Up",       "Up",        "Up",
  "ATP7A",   "Copper transport; brain metabolism (reduced in MDD)", "Up",       "Down",      "Down",
  "CTNND2",  "Neurodevelopmental risk; cell adhesion",              "Up",       "Up",        "Up",
  "CHRNA7",  "Cognitive symptoms; cholinergic signaling",           "Up",       "Up",        "Up",
  "SEMA3A",  "Axon guidance; synaptic plasticity in MDD",           "Up",       "Up",        "Up"
)
write_csv(fig7a_tbl, "Figure_7A_gene_table.csv")

# Render the table as a figure panel (gridExtra/ggpubr if available)
if (requireNamespace("ggpubr", quietly = TRUE)) {
  arrow_col <- function(x) ifelse(x == "Up", "#B2182B",
                            ifelse(x == "Down", "#2166AC", "grey50"))
  thm <- ggpubr::ttheme(
    base_size = 9,
    colnames.style = ggpubr::colnames_style(fill = "#0072B2", color = "white"),
    tbody.style    = ggpubr::tbody_style(fill = c("white", "#F2F2F2")))
  p7a <- ggpubr::ggtexttable(fig7a_tbl, rows = NULL, theme = thm)
  ggplot2::ggsave("Figure_7A_table.tiff", p7a, width = 9, height = 3.2,
                  dpi = 300, compression = "lzw")
} else {
  message("[7A] install 'ggpubr' to render the gene table as a TIFF; ",
          "table data written to Figure_7A_gene_table.csv.")
}

## ---- 7B. hsa-miR-513a-3p effect across timepoints (DESeq2 estimates) ------ ##
if (file.exists("ColData_Final.csv") && file.exists("MDD_Final_Counts.csv") &&
    requireNamespace("DESeq2", quietly = TRUE)) {

  suppressPackageStartupMessages(library(DESeq2))

  metadata   <- read.csv("ColData_Final.csv", fileEncoding = "UTF-8-BOM")
  Count_Data <- read.csv("MDD_Final_Counts.csv", row.names = 1,
                         header = TRUE, check.names = FALSE)
  rownames(metadata) <- metadata$Sample_ID
  metadata <- metadata[colnames(Count_Data), , drop = FALSE]
  metadata$ID        <- factor(metadata$ID)
  metadata$Treatment <- factor(metadata$Treatment,
                               levels = c("PRE","MID","POST","75 min POST"))

  dds <- DESeqDataSetFromMatrix(Count_Data, metadata, ~ ID + Treatment)
  dds <- estimateSizeFactors(dds)
  dds <- dds[rowSums(counts(dds, normalized = TRUE) >= 10) >= 3, ]
  dds$Treatment <- relevel(dds$Treatment, ref = "PRE")
  dds <- DESeq(dds)
  ddsClean <- dds[which(mcols(dds)$betaConv), ]

  g <- "hsa-miR-513a-3p"
  grab <- function(tp) {
    r <- results(ddsClean, contrast = c("Treatment", tp, "PRE"))[g, ]
    data.frame(Timepoint = tp, log2FC = r$log2FoldChange,
               lfcSE = r$lfcSE, padj = r$padj)
  }
  est <- rbind(grab("MID"), grab("POST"), grab("75 min POST"))
  est$Timepoint <- factor(c("Mid-AAE","Post-AAE","75 min Post-AAE"),
                          levels = c("Mid-AAE","Post-AAE","75 min Post-AAE"))
  est$lab <- sprintf("log2FC = %.2f\nFDR = %.3f", est$log2FC, est$padj)

  fc_fill <- c("Mid-AAE" = unname(tp_colors["MID"]),
               "Post-AAE" = unname(tp_colors["POST"]),
               "75 min Post-AAE" = unname(tp_colors["75 min POST"]))

  fig7b <- ggplot(est, aes(Timepoint, log2FC, fill = Timepoint)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(aes(ymin = log2FC - lfcSE, ymax = log2FC + lfcSE),
                  width = 0.18, linewidth = 0.5) +
    geom_text(aes(label = lab, y = log2FC + lfcSE + 0.2),
              size = 3, fontface = "italic", vjust = 0) +
    scale_fill_manual(values = fc_fill, guide = "none") +
    labs(title = "hsa-miR-513a-3p Differential Expression vs Pre-AAE",
         subtitle = "DESeq2 model estimate (~ ID + Treatment); error bars = lfcSE",
         x = NULL, y = expression(log[2]~"Fold Change vs Pre-AAE")) +
    theme_pub + theme(plot.subtitle = element_text(size = 9))

  ggsave("Figure_7B_513a_log2FC.tiff", fig7b, width = 6.5, height = 5,
         dpi = 300, compression = "lzw")
} else {
  message("[7B] needs ColData_Final.csv, MDD_Final_Counts.csv and DESeq2 - skipping.")
}

## ---- 7C. Bubble plot: top-10 enriched pathways for hsa-miR-513a-3p -------- ##
top10_513a <- data.frame(
  GO_Term = c(
    "synapse assembly",
    "regulation of cellular response to transforming growth factor beta stimulus",
    "negative regulation of TOR signaling",
    "regulation of transforming growth factor beta receptor signaling pathway",
    "cellular response to leukemia inhibitory factor",
    "regulation of translational initiation",
    "negative regulation of mRNA metabolic process",
    "response to leukemia inhibitory factor",
    "miRNA processing",
    "dendrite morphogenesis"),
  FoldEnrichment = c(2.44, 2.21, 2.69, 2.15, 2.31, 2.42, 2.29, 2.27, 3.39, 2.62),
  Bonferroni     = c(0.0000571, 0.00188, 0.00307, 0.00566, 0.0186,
                     0.0219, 0.0229, 0.0284, 0.0300, 0.0358))
top10_513a$GO_Term <- factor(top10_513a$GO_Term, levels = rev(top10_513a$GO_Term))

fig7c <- ggplot(top10_513a, aes(x = "", y = GO_Term)) +
  geom_point(aes(size = FoldEnrichment, color = -log10(Bonferroni))) +
  scale_color_viridis_c(option = viridis_option, name = expression(-log[10]~"(Bonferroni)")) +
  scale_size_continuous(name = "Fold Enrichment", range = c(3, 9)) +
  labs(title = "Top 10 Enriched Pathways for hsa-miR-513a-3p Targets",
       x = NULL, y = "GO Biological Process") +
  theme_pub +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggsave("Figure_7C_513a_bubble.tiff", fig7c, width = 8, height = 5,
       dpi = 300, compression = "lzw")

message("\nDone. Figure 6 & 7 TIFFs written to: ", getwd())
