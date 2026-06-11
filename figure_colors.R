
## ----- Timepoint colors (Okabe-Ito) ------------------------------------- ##
## Pre  = black, Mid = bluish-green, Post = blue, 75 min = vermillion/orange.
tp_colors <- c(
  "PRE"         = "#000000",  # black
  "MID"         = "#009E73",  # bluish green
  "POST"        = "#0072B2",  # blue
  "75 min POST" = "#D55E00"   # vermillion
)

## Pretty labels (for legends / facets)
tp_labels <- c(
  "PRE"         = "Pre-AAE",
  "MID"         = "Mid-AAE",
  "POST"        = "Post-AAE",
  "75 min POST" = "75 minutes post-AAE"
)

## Same colors keyed by the short pattern labels used in annotations
tp_colors_short <- c(
  "MID"       = "#009E73",
  "POST"      = "#0072B2",
  "75MIN"     = "#D55E00",
  "75minPOST" = "#D55E00"
)

## ----- DEG-pattern annotation colors (Fig 4A row annotation) ------------- ##
## Singles reuse the timepoint colors; combos use blended Okabe-Ito tones.
deg_pattern_colors <- c(
  "Mid-AAE Only"     = "#009E73",
  "Post-AAE Only"    = "#0072B2",
  "75min-AAE Only"   = "#D55E00",
  "Mid & Post-AAE"   = "#56B4E9",  # sky blue
  "Post & 75min-AAE" = "#CC79A7",  # reddish purple
  "Mid & 75min-AAE"  = "#E69F00"   # orange
)

## ----- Continuous (bubble plots, -log10 p) ------------------------------- ##
## Use viridis. In ggplot:  scale_color_viridis_c(option = "viridis")
viridis_option <- "viridis"

## ----- Diverging palette for log2FC heatmaps (4A, 4B, 7b) ---------------- ##
## One consistent scale across every heatmap. Symmetric around 0.
divergent_palette <- function(n = 100) {
  colorRampPalette(c("#2166AC",  # strong blue  (down)
                     "#67A9CF",
                     "#F7F7F7",  # white        (no change)
                     "#EF8A62",
                     "#B2182B"   # strong red   (up)
                     ))(n)
}
## Shared symmetric breaks helper (cap controls the saturation range)
divergent_breaks <- function(cap = 10, n = 100) seq(-cap, cap, length.out = n + 1)

## ----- Sequential palette for similarity heatmaps (Fig 5B Jaccard) ------- ##
jaccard_palette <- function(n = 100)
  colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(n)  # white -> blue (CB-safe)

message("Loaded shared color-blind-safe palette (figure_colors.R).")
