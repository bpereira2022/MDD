# General QC Analysis

## ===========================================================================
## Setting Up (:
## ===========================================================================
# Install 
if (!requireNamespace("DESeq2")) install.packages("BiocManager"); BiocManager::install("DESeq2")
if (!requireNamespace("pheatmap")) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer")) install.packages("RColorBrewer")
if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("vsn")) BiocManager::install("vsn")

# Load packages
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(vsn)

## ===========================================================================
## General QC
## ===========================================================================

# Create DESeqDataSet

metadata <- read.csv("ColData_MDD.csv")
Count_Data <- read.csv("MDD_Final_Counts.csv", row.names = 1, header = TRUE, check.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = Count_Data,
                              colData = metadata,
                              design = ~ 1)  # no model yet, just QC

# Filter out low count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# Library size check
colSums(counts(dds))
barplot(colSums(counts(dds)), las=2, main="Library sizes", ylab="Total counts")

# Sample clustering and PCA
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)  # variance-stabilizing transformation

# Flor Plotting: Perform VST transformation
vst_counts <- assay(varianceStabilizingTransformation(dds))

# Calculate library sizes (sum of VST values per sample)
vst_lib_sizes <- colSums(vst_counts)

# Plot
barplot(vst_lib_sizes, las=2, main="VST-Transformed: Library Sizes", ylab="Sum of VST values")

# Calculate library sizes (total counts per sample)
raw_lib_sizes <- colSums(counts(dds))

# Combine into a matrix for side-by-side plotting
library_sizes <- rbind(raw_lib_sizes, vst_lib_sizes)
rownames(library_sizes) <- c("Raw", "VST")

barplot(library_sizes, beside=TRUE, las=2,
        main="Library Sizes: Raw vs VST",
        ylab="Library size / Sum of VST values",
        col=c("skyblue", "orange"))
legend("topright", legend=c("Raw", "VST"), fill=c("skyblue", "orange"))


# PCA
pcaData <- plotPCA(vsd, intgroup=colnames(metadata), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x=PC1, y=PC2, color=metadata[[1]])) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds)
colnames(sampleDistMatrix) <- colnames(dds)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255))

# Mean-variance relationship & Negative Binomial fit
dds_fit <- estimateSizeFactors(dds)
dds_fit <- estimateDispersions(dds_fit)
plotDispEsts(dds_fit)  # dispersion plot
# Red line: expected dispersion under the negative binomial
# Black points: gene-wise estimates
# Blue: final dispersion used

# Mean-sd plot (diagnosing heteroscedasticity) ===
meanSdPlot(assay(vsd))  # from vsn package

# Check for outliers using Cook's distance
dds_fit <- nbinomWaldTest(dds_fit)  # doesn't do testing, just needed to get mcols
cooks <- assays(dds_fit)[["cooks"]]
boxplot(log10(cooks + 1), las=2, main="Cook's distance (log10)")

# Top expressed genes
# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Use variance-stabilized data
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Get the normalized expression matrix
expr_matrix <- assay(vsd)

# Calculate mean expression for each miRNA
mean_expression <- rowMeans(expr_matrix)

# Create a data frame with miRNA names and mean expression
top_miRNAs_df <- data.frame(
  miRNA = rownames(expr_matrix),
  mean_expr = mean_expression
)

# Get top 20 most expressed miRNAs
top20 <- top_miRNAs_df %>%
  arrange(desc(mean_expr)) %>%
  slice_head(n = 20)

# Reorder
top20$miRNA <- factor(top20$miRNA, levels = top20$miRNA[order(top20$mean_expr)])

# Plot the top 20
ggplot(top20, aes(x = miRNA, y = mean_expr)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Most Expressed miRNAs",
       x = "miRNA",
       y = "Mean Expression (VST normalized)") +
  theme_minimal(base_size = 14)

## ===========================================================================
## Outlier Detection
## ===========================================================================
# Using isotree
library(isotree)

dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3

# Running our model
dds <- dds[idx,]
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
sample_matrix <- t(norm_counts)

# Apply Isolation Forest
iso_model <- isolation.forest(sample_matrix, ntrees = 100)
outlier_scores <- predict(iso_model, sample_matrix, type = "score")

# Label outliers (e.g., top 5%)
threshold <- quantile(outlier_scores, 0.95)
outlier_samples <- which(outlier_scores >= threshold)

# View results
outlier_sample_names <- rownames(sample_matrix)[outlier_samples]

# Print names
outlier_sample_names

# got "352_E1" as > 0.5

# Cooks Distance
dds_fit <- estimateSizeFactors(dds)
dds_fit <- estimateDispersions(dds_fit)
# Check for outliers using Cook's distance
dds_fit <- nbinomWaldTest(dds_fit)  # doesn't do testing, just needed to get mcols
cooks <- assays(dds_fit)[["cooks"]]
boxplot(log10(cooks + 1), las=2, main="Cook's distance (log10)")

# PCA Plot
# Sample clustering and PCA
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)  # variance-stabilizing transformation

# Attach sample metadata to PCA data
pcaData <- plotPCA(vsd, intgroup = colnames(metadata), returnData = TRUE)

# Add color group column from metadata (e.g., group by Stress, Run, etc.)
pcaData$Group <- metadata[[1]]  # or change this to metadata$Run or metadata$Stress

# Calculate variance from origin in PCA space
pcaData$PC_variance <- with(pcaData, PC1^2 + PC2^2)

# Select top 5 most variant samples
top10_samples <- top10_samples <- pcaData[order(-pcaData$PC_variance), ][1:10, ]

# Optional: add sample names if rownames are not informative
pcaData$Sample <- rownames(pcaData)
top10_samples$Sample <- rownames(top10_samples)

library(ggplot2)
library(ggrepel)

ggplot(pcaData, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(
    data = top10_samples,
    aes(label = Sample),  # use Sample column, not rownames(top5_samples)
    size = 3.5,
    box.padding = 0.4
  ) +
  xlab(paste0("PC1: ", round(100 * attr(pcaData, "percentVar")[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaData, "percentVar")[2]), "% variance")) +
  theme_bw()

# Cooks & PCA did not confirm 352_E1 as an global outlier, so all samples were kept for analysis.
