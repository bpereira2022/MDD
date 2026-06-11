
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(viridis); library(stringr)
})
select <- dplyr::select; filter <- dplyr::filter

dir.create("results", showWarnings = FALSE)
out <- function(f) file.path("results", f)
source("figure_colors.R")   # for viridis_option and category palette consistency

XLSX  <- "PANTHER_PATHWAY_RESULTS.xlsx"   # put in working dir
TOP_N <- 10                                # terms shown per timepoint in 6A

## ---- 1. Read all three sheets into one tidy frame ------------------------- ##
read_sheet <- function(sheet, label) {
  d <- read_excel(XLSX, sheet = sheet)
  # PANTHER columns: term | #ref | #list | expected | Fold Enrichment | +/- | P value
  names(d)[1] <- "Term"
  fe  <- grep("fold", names(d), ignore.case = TRUE, value = TRUE)[1]
  pv  <- grep("^p.?value|FDR|adjust", names(d), ignore.case = TRUE, value = TRUE)[1]
  # the genes-in-list count is the 3rd column (second '#')
  gcol <- names(d)[3]
  tibble(Timepoint = label,
         Term      = d$Term,
         genes     = as.numeric(d[[gcol]]),
         FoldEnrichment = as.numeric(d[[fe]]),
         Pvalue    = as.numeric(d[[pv]])) %>%
    filter(!is.na(Term), !is.na(FoldEnrichment))
}
sheets <- readxl::excel_sheets(XLSX)
lab_map <- c("MID" = "Mid-AAE", "POST" = "Post-AAE", "75 min POST" = "75 min Post-AAE")
enr <- bind_rows(lapply(sheets, function(s)
  read_sheet(s, ifelse(s %in% names(lab_map), lab_map[[s]], s))))
enr$Timepoint <- factor(enr$Timepoint,
                        levels = c("Mid-AAE","Post-AAE","75 min Post-AAE"))
## CHECKPOINT 1
cat("Loaded PANTHER results:\n")
print(enr %>% group_by(Timepoint) %>%
        summarise(terms = n(), FE_min = round(min(FoldEnrichment),2),
                  FE_max = round(max(FoldEnrichment),2),
                  min_genes = min(genes)))

## ---- 2. Umbrella categories (keyword rules, ORDER = priority) ------------- ##
## Specific neuronal-function terms first, then structural, then RNA, etc.
## First match wins. Grounded in the actual top terms; edit as you see fit.
categorize <- function(term) {
  t <- tolower(term)
  dplyr::case_when(
    str_detect(t, "nerve impulse|synaptic transmission|trans-synaptic|neurotransmitter|action potential|membrane potential") ~ "Synaptic Transmission",
    str_detect(t, "dendrit|spine|postsynap|\\bsynapse\\b|synaptic|filopodium") ~ "Dendrite & Synapse",
    str_detect(t, "axon|neuron projection|projection neuron|fascicul")        ~ "Axon & Projection",
    str_detect(t, "learning|memory|cognition|behavio")                        ~ "Learning & Behavior",
    str_detect(t, "migrat|neurogenesis|gliogenesis|forebrain|telencephalon|cerebral|cortex|neuron|neural|nervous system|brain") ~ "Neuronal Development",
    str_detect(t, "\\brna\\b|mrna|splic|translat|catabolic|nucleobase|nucleotide|chromatin") ~ "RNA & Chromatin",
    str_detect(t, "ion transport|ion transmembrane|potassium|sodium|calcium|cation|metal ion|transmembrane transport") ~ "Ion Transport",
    str_detect(t, "adhesion|cell junction|cell-cell|polarity")                ~ "Adhesion & Junction",
    str_detect(t, "signal|kinase|hippo|tgf|transforming growth|hormone")      ~ "Cell Signaling",
    str_detect(t, "develop|morphogen|differentiat|growth|angiogen|vascul")    ~ "Development & Growth",
    # broad top-level GO nodes (roots of the ontology) - named rather than dumped
    str_detect(t, "regulation of (biological|cellular|macromolecule|primary metabolic|metabolic|response)|metabolic process|cellular component organization|multicellular organismal process|^biological|^cellular process|biological_process") ~ "General Regulation",
    TRUE ~ "Other"
  )
}
enr <- enr %>% mutate(Category = categorize(Term))
write_csv(enr, out("panther_categorized.csv"))
## CHECKPOINT 2: category assignment overview (full list)
cat("\nCategory assignment (ALL significant terms):\n")
print(enr %>% count(Category, sort = TRUE))
cat("Note: the full list contains many generic regulatory terms; Figure 6B is\n",
    "built from the TOP terms per timepoint (below), where 'Other' is minimal.\n")

## ---- 3. Figure 6A: top-N terms per timepoint by fold enrichment ----------- ##
top <- enr %>% group_by(Timepoint) %>%
  arrange(desc(FoldEnrichment), .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>% ungroup()
write_csv(top, out("bubble_top_all.csv"))

# order y-axis by fold enrichment within timepoint
term_order <- top %>% arrange(Timepoint, FoldEnrichment) %>% pull(Term) %>% unique()
top$Term <- factor(top$Term, levels = term_order)

fig6a <- ggplot(top, aes(Timepoint, Term)) +
  geom_point(aes(size = FoldEnrichment, color = -log10(Pvalue))) +
  scale_color_viridis_c(option = viridis_option, name = expression(-log[10]~"(P value)")) +
  scale_size_continuous(name = "Fold Enrichment", range = c(2.5, 8)) +
  labs(title = "Top Enriched GO Biological Processes by Timepoint",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank())
ggsave(out("Figure_6A_bubble.tiff"), fig6a, width = 9, height = 8,
       dpi = 300, compression = "lzw")
cat("\nWrote", out("Figure_6A_bubble.tiff"), "\n")

## ---- 4. Figure 6B: umbrella-category composition per timepoint ------------ ##
N_CAT <- 20   # number of top terms (by fold enrichment) used for composition
top_cat <- enr %>% group_by(Timepoint) %>%
  arrange(desc(FoldEnrichment), .by_group = TRUE) %>%
  slice_head(n = N_CAT) %>% ungroup()
cat_counts <- top_cat %>% count(Timepoint, Category)
write_csv(cat_counts, out("category_counts.csv"))
## CHECKPOINT 3: confirm 'Other' is small in the top-N composition
other_top <- top_cat %>% group_by(Timepoint) %>%
  summarise(pct_other = round(100 * mean(Category == "Other"), 0))
cat("\n'Other' share within top", N_CAT, "terms per timepoint:\n"); print(other_top)

cat_levels <- c("RNA & Chromatin","Dendrite & Synapse","Axon & Projection",
                "Synaptic Transmission","Neuronal Development","Learning & Behavior",
                "Ion Transport","Adhesion & Junction","Cell Signaling",
                "Development & Growth","General Regulation","Other")
cat_colors <- c("RNA & Chromatin"        = "#0072B2",  # blue
                "Dendrite & Synapse"     = "#009E73",  # green
                "Axon & Projection"      = "#56B4E9",  # sky blue
                "Synaptic Transmission"  = "#CC79A7",  # reddish purple
                "Neuronal Development"   = "#E69F00",  # orange
                "Learning & Behavior"    = "#F0E442",  # yellow
                "Ion Transport"          = "#999999",  # grey
                "Adhesion & Junction"    = "#661100",  # dark red
                "Cell Signaling"         = "#D55E00",  # vermillion
                "Development & Growth"   = "#117733",  # dark green
                "General Regulation"     = "#AA4499")  # purple
cat_counts$Category <- factor(cat_counts$Category, levels = cat_levels)

fig6b <- ggplot(cat_counts, aes(Timepoint, n, fill = Category)) +
  geom_col(position = "fill", width = 0.7) +     # proportional composition
  scale_fill_manual(values = cat_colors, name = "Pathway Category", drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = paste0("Enriched-Pathway Composition by Timepoint (top ",
                      N_CAT, " terms)"),
       x = NULL, y = "Proportion of top GO terms") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(out("Figure_6B_categories.tiff"), fig6b, width = 8, height = 5,
       dpi = 300, compression = "lzw")
cat("Wrote", out("Figure_6B_categories.tiff"), "\n")

cat("\nDONE. Review panther_categorized.csv to check category assignments;\n",
    "adjust the keyword rules in categorize() if any terms are miscategorized.\n")

## ---- 4. Figure 7B: miR-1252-3p Bubble Plot ------------ ##
library(ggplot2)
library(viridis)
library(dplyr)

# Read your pathway results
pathways <- read.csv("special_miR_pathways.csv")

# Convert Pvalue to numeric (in case Excel saved it in scientific notation)
pathways$Pvalue <- as.numeric(pathways$Pvalue)

# Sort by Fold Enrichment and keep top pathways
top <- pathways %>%
  arrange(desc(FoldEnrichment))

# Order terms by fold enrichment
term_order <- top %>% 
  arrange(FoldEnrichment) %>% 
  pull(Term) %>% 
  unique()

top$Term <- factor(top$Term, levels = term_order)

# Create bubble plot
fig_pathways <- ggplot(top, aes(x = 1, y = Term)) +
  geom_point(aes(size = FoldEnrichment, color = -log10(Pvalue))) +
  scale_color_viridis_c(option = "viridis", name = expression(-log[10]~"(P value)")) +
  scale_size_continuous(name = "Fold Enrichment", range = c(3, 10)) +
  labs(title = "Enriched GO Biological Processes",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave("pathway_enrichment_bubble.pdf", fig_pathways, width = 7, height = 4, dpi = 300)
ggsave("pathway_enrichment_bubble.tiff", fig_pathways, width = 7, height = 4, dpi = 300, compression = "lzw")

cat("\n✓ Bubble plot saved!\n")
print(fig_pathways)