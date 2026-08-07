library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

# Run after getting dds object

norm_counts <- counts(ddsClean, normalized = TRUE)

metadata_qc <- as.data.frame(colData(ddsClean)) %>%
  rownames_to_column("SampleID") %>%
  rename(SubjectID = ID, Timepoint = Treatment)

timepoint_order <- c("PRE", "MID", "POST", "75 min POST")

hemolysis_numerator   <- c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p")
hemolysis_denominator <- c("hsa-miR-23a-3p", "hsa-miR-93-5p")

platelet_mirnas <- c(
  "hsa-miR-223-3p",
  "hsa-miR-126-3p",
  "hsa-miR-21-5p",
  "hsa-miR-92a-1-5p",
  "hsa-miR-191-5p"
)

# =========================================================================
# HEMOLYSIS INDEX: computation, Friedman test, boxplot, and PCA check
# =========================================================================

num_present <- hemolysis_numerator[hemolysis_numerator %in% rownames(norm_counts)]
den_present <- hemolysis_denominator[hemolysis_denominator %in% rownames(norm_counts)]

cat("Hemolysis numerator miRNAs found:", paste(num_present, collapse = ", "), "\n")
cat("Hemolysis denominator miRNAs found:", paste(den_present, collapse = ", "), "\n")

primary_num <- if ("hsa-miR-451a"   %in% num_present) "hsa-miR-451a"   else num_present[1]
primary_den <- if ("hsa-miR-23a-3p" %in% den_present) "hsa-miR-23a-3p" else den_present[1]

cat("Using ratio:", primary_num, "/", primary_den, "\n")

numerator_counts   <- as.numeric(norm_counts[primary_num, ]) + 1
denominator_counts <- as.numeric(norm_counts[primary_den, ]) + 1
hemolysis_index    <- log2(numerator_counts / denominator_counts)

hi_df <- data.frame(
  SampleID       = colnames(norm_counts),
  HemolysisIndex = hemolysis_index
) %>%
  left_join(metadata_qc, by = "SampleID") %>%
  mutate(Timepoint = factor(Timepoint, levels = timepoint_order))

hi_avg <- hi_df %>%
  group_by(SubjectID, Timepoint) %>%
  summarise(HemolysisIndex = mean(HemolysisIndex, na.rm = TRUE), .groups = "drop")

complete_subjects <- hi_avg %>%
  group_by(SubjectID) %>%
  filter(n() == length(timepoint_order)) %>%
  ungroup()

complete_subjects$Timepoint <- droplevels(complete_subjects$Timepoint)
complete_subjects$SubjectID <- droplevels(complete_subjects$SubjectID)

cat("Subjects with complete hemolysis data across all timepoints:",
    n_distinct(complete_subjects$SubjectID), "\n")

friedman_result <- friedman.test(
  HemolysisIndex ~ Timepoint | SubjectID,
  data = complete_subjects
)

cat("\nFriedman test for hemolysis index across timepoints:\n")
print(friedman_result)

if (friedman_result$p.value < 0.05) {
  cat("\nWARNING: Hemolysis index differs significantly across timepoints.\n")
} else {
  cat("\nHemolysis index does not differ significantly across timepoints (p =",
      round(friedman_result$p.value, 3), ").\n")
  cat("Gross hemolysis is unlikely to account for observed temporal miRNA patterns.\n")
}

hemolysis_threshold <- log2(7)

p_hemolysis <- ggplot(hi_df, aes(x = Timepoint, y = HemolysisIndex)) +
  geom_hline(yintercept = hemolysis_threshold,
             linetype = "dashed", color = "#CC3333", linewidth = 0.8) +
  annotate("text", x = "PRE", y = hemolysis_threshold + 0.15,
           label = paste0("Hemolysis threshold (log2 = ",
                          round(hemolysis_threshold, 2), ")"),
           size = 3, color = "#CC3333", hjust = 0) +
  geom_boxplot(aes(fill = Timepoint), width = 0.5, alpha = 0.7,
               outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "grey30") +
  scale_fill_manual(values = c(
    "PRE"         = "#AAAAAA",
    "MID"         = "#009E73",
    "POST"        = "#0072B2",
    "75 min POST" = "#D55E00"
  )) +
  labs(
    title    = "Hemolysis Index Across Timepoints",
    subtitle = paste0("log2(", primary_num, " / ", primary_den, ")"),
    x        = "Timepoint",
    y        = "Hemolysis Index (log2 ratio)",
    caption  = paste0("Friedman test p = ", round(friedman_result$p.value, 3),
                      " | Dashed line = hemolysis threshold")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
    plot.caption    = element_text(color = "grey50", size = 9)
  )

# ----- PCA check colored by hemolysis index -----

vsd <- vst(ddsClean, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "Treatment", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_data <- pca_data %>%
  rownames_to_column("SampleID") %>%
  left_join(hi_df %>% select(SampleID, HemolysisIndex), by = "SampleID")

p_pca_hemolysis <- ggplot(pca_data, aes(x = PC1, y = PC2, color = HemolysisIndex)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Hemolysis\nIndex") +
  labs(
    title = "PCA of miRNA Expression Colored by Hemolysis Index",
    x     = paste0("PC1: ", percentVar[1], "% variance"),
    y     = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

hemolysis_combined <- p_hemolysis / p_pca_hemolysis +
  plot_annotation(
    title    = "Post-hoc Hemolysis Assessment",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave(out("Supplementary_Figure2_Hemolysis_QC.tiff"),
       plot = hemolysis_combined, width = 8, height = 12, dpi = 300,
       compression = "lzw")

ggsave(out("Supplementary_Figure2_Hemolysis_QC.pdf"),
       plot = hemolysis_combined, width = 8, height = 12)

write.csv(
  hi_df %>% select(SampleID, Timepoint, SubjectID, HemolysisIndex),
  out("Supplementary_Table_HemolysisIndex.csv"),
  row.names = FALSE
)

# =========================================================================
# PLATELET-ENRICHED miRNAs: per-miRNA Friedman tests and boxplot panel
# =========================================================================

plt_present <- platelet_mirnas[platelet_mirnas %in% rownames(norm_counts)]
cat("\nPlatelet miRNAs found in data:", paste(plt_present, collapse = ", "), "\n")

plt_long <- as.data.frame(norm_counts[plt_present, , drop = FALSE]) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(-miRNA, names_to = "SampleID", values_to = "NormCount") %>%
  left_join(metadata_qc, by = "SampleID") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoint_order),
    DE_flag   = ifelse(miRNA %in% c("hsa-miR-21-5p", "hsa-miR-92a-1-5p"),
                       "Also DE in study", "Not DE")
  )

# Friedman test for each platelet-enriched miRNA across timepoints

plt_avg <- plt_long %>%
  group_by(miRNA, SubjectID, Timepoint) %>%
  summarise(NormCount = mean(NormCount, na.rm = TRUE), .groups = "drop")

platelet_friedman_results <- lapply(plt_present, function(mir) {

  df_mir <- plt_avg %>%
    filter(miRNA == mir) %>%
    group_by(SubjectID) %>%
    filter(n() == length(timepoint_order)) %>%
    ungroup()

  df_mir$Timepoint <- droplevels(df_mir$Timepoint)
  df_mir$SubjectID  <- droplevels(df_mir$SubjectID)

  n_subj <- n_distinct(df_mir$SubjectID)

  if (n_subj < 3 || nlevels(df_mir$Timepoint) < length(timepoint_order)) {
    return(data.frame(
      miRNA      = mir,
      n_subjects = n_subj,
      statistic  = NA_real_,
      df         = NA_real_,
      p_value    = NA_real_,
      note       = "Insufficient complete cases across all timepoints"
    ))
  }

  ft <- friedman.test(NormCount ~ Timepoint | SubjectID, data = df_mir)

  data.frame(
    miRNA      = mir,
    n_subjects = n_subj,
    statistic  = unname(ft$statistic),
    df         = unname(ft$parameter),
    p_value    = ft$p.value,
    note       = NA_character_
  )
}) %>% bind_rows()

platelet_friedman_results <- platelet_friedman_results %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cat("\nFriedman test results for platelet-enriched miRNAs across timepoints:\n")
print(platelet_friedman_results)

write.csv(
  platelet_friedman_results,
  out("Supplementary_Table_PlateletMiRNA_FriedmanTests.csv"),
  row.names = FALSE
)

p_platelet <- ggplot(plt_long, aes(x = Timepoint, y = log2(NormCount + 1))) +
  geom_boxplot(aes(fill = Timepoint), width = 0.5, alpha = 0.7,
               outlier.shape = NA) +
  geom_jitter(aes(color = DE_flag), width = 0.1, size = 1.5, alpha = 0.7) +
  facet_wrap(~ miRNA, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "PRE"         = "#AAAAAA",
    "MID"         = "#009E73",
    "POST"        = "#0072B2",
    "75 min POST" = "#D55E00"
  )) +
  scale_color_manual(
    values = c("Also DE in study" = "#CC3333", "Not DE" = "#333333"),
    name   = ""
  ) +
  labs(
    title    = "Platelet-Enriched miRNAs Across Timepoints",
    subtitle = "Normalized counts (log2 + 1) | Red = also differentially expressed in main analysis",
    x        = "Timepoint",
    y        = "log2(Normalized Count + 1)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 9),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  ) +
  plot_annotation(
    title    = "Post-hoc Platelet Contamination Assessment",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave(out("Supplementary_Figure3_Platelet_QC.tiff"),
       plot = p_platelet, width = 10, height = 8, dpi = 300,
       compression = "lzw")

ggsave(out("Supplementary_Figure3_Platelet_QC.pdf"),
       plot = p_platelet, width = 10, height = 8)

cat("\nDone. Outputs saved to results/ folder.\n")
