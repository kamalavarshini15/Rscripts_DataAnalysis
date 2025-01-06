library(DESeq2)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(ggrepel)
setwd("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq")
# Load count data from CSV files
control_counts <- read_csv("COPD_control_count.csv")
copd_counts <- read_csv("COPD_case_count.csv")
ra_counts <- read_csv("RA_count.csv")

# Merge count data (assuming first column is 'Gene' in each file)
merged_counts <- reduce(list(control_counts, copd_counts, ra_counts), full_join, by = "Gene")
head(merged_counts)
sample_names <- colnames(merged_counts)

gene_names <- merged_counts$Gene
merged_counts <- merged_counts %>% select(-Gene) 
# set gene column to rownames
merged_counts <- merged_counts %>%column_to_rownames(var = "Gene")
colnames(merged_counts) <- gsub("_gene_counts", "", colnames(merged_counts))
write.csv(merged_counts,file="Deseq_dataset.csv")
# Load sample information
sample_info <- read_csv("sample_info.csv")  # Load your sample information CSV
rownames(sample_info) <- sample_info$Sample

# Ensure columns in counts match samples in sample_info
counts <- merged_counts[,rownames(sample_info)]

# Create DESeq2 dataset using sample condition from sample_info
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ Condition)
# View the first 10 row names to check for gene names
head(dds)
# Run the DESeq2 pipeline
dds <- DESeq(dds)
colnames(dds)[1:10]
# Extract normalized counts (optional step, DESeq2 internally uses normalized data for analysis)
normalized_counts <- counts(dds, normalized = TRUE)

# Print normalized counts (optional, just to see how they look)
head(normalized_counts)







# Function to create volcano plot
plot_volcano <- function(res, title , gene_names) {
  res$logp <- -log10(res$pvalue)
  res$significance <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 2, "Significant", "Not Significant")
  ggplot(res, aes(x = log2FoldChange, y = logp, color = significance)) +
    geom_point(alpha = 0.8) + scale_color_manual(values = c("grey", "red")) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(p-value)") + theme_minimal()
}


plot_volcano <- function(res, title) {
  res$logp <- -log10(res$pvalue)
  res$significance <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0, "Upregulated",
                             ifelse(res$padj < 0.05 & res$log2FoldChange < 0, "Downregulated", "Not Significant"))
  
  # Create a volcano plot
  p <- ggplot(res, aes(x = log2FoldChange, y = logp, color = significance)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal()
  
  # Add gene names to significant points
  p + geom_text_repel(data = subset(res, significance != "Not Significant"),
                      aes(label = gene), 
                      size = 3, 
                      max.overlaps = Inf)  # Avoid overlapping labels
}

# Get results and plot for Control vs COPD
res_copd <- results(dds, contrast = c("Condition", "COPD", "Control"))
res_copd_df <- as.data.frame(res_copd)

# Add gene names as a new column
res_copd_df$gene <- gene_names

plot_volcano(res_copd, "COPD vs Control")

# Get results and plot for Control vs RA
res_ra <- results(dds, contrast = c("Condition", "RA", "Control"))


# Convert the results to a data frame
res_ra_df <- as.data.frame(res_ra)

# Add gene names as a new column
res_ra_df$gene <- gene_names
plot_volcano(res_ra, "RA vs Control")

# Get results and plot for COPD vs RA
res_copd_ra <- results(dds, contrast = c("Condition", "RA", "COPD"))

res_copd_ra_df <- as.data.frame(res_copd_ra)

# Add gene names as a new column
res_copd_ra_df$gene <- gene_names
plot_volcano(res_copd_ra, "RA vs COPD")

plot_volcano <- function(res, title) {
  res$logp <- -log10(res$pvalue)
  res$significance <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0, "Upregulated",
                             ifelse(res$padj < 0.05 & res$log2FoldChange < 0, "Downregulated", "Not Significant"))
  
  # Create a volcano plot with color coding for upregulated and downregulated genes
  p <- ggplot(res, aes(x = log2FoldChange, y = logp, color = significance)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal()
  
  # Add gene names to significant points
  #p + geom_text_repel(data = subset(res, significance != "Not Significant"),
                      #aes(label = gene), 
                      #size = 3, 
                      #max.overlaps = Inf)  # Avoid overlapping labels
}

# Get results and plot for Control vs COPD
res_copd <- results(dds, contrast = c("Condition", "COPD", "Control"))
res_copd_df <- as.data.frame(res_copd)

# Add gene names as a new column
res_copd_df$gene <- gene_names
print(res_copd_df)

plot_volcano(res_copd_df, "COPD vs Control")

# Get results and plot for Control vs RA
res_ra <- results(dds, contrast = c("Condition", "RA", "Control"))


# Convert the results to a data frame
res_ra_df <- as.data.frame(res_ra)

# Add gene names as a new column
res_ra_df$gene <- gene_names
plot_volcano(res_ra_df, "RA vs Control")













library(ggrepel)  # Ensure this library is loaded for text repelling

plot_volcano <- function(res, title) {
  res$logp <- -log10(res$pvalue)
  res$significance <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")
  
  # Create a volcano plot
  p <- ggplot(res, aes(x = log2FoldChange, y = logp, color = significance)) +
    geom_point(alpha = 0.8) + 
    scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
    #scale_color_manual(values = c("grey", "red")) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal()
  
  # Add gene names to significant points
  p + geom_text_repel(data = subset(res, significance == "Significant"),
                      aes(label = gene), 
                      size = 3, 
                      max.overlaps = Inf)  # Avoid overlapping labels
}

res_copd <- results(dds, contrast = c("Condition", "COPD", "Control"))
head(res_copd)
res_copd_df <- as.data.frame(res_copd)

# Add gene names as a new column
res_copd_df$gene <- gene_names




# Plot the volcano plot with gene names
plot_volcano(res_copd_df, "COPD vs Control")
significant_genes <- res_copd_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  select(gene, log2FoldChange, padj)
print(significant_genes)
write.csv(significant_genes, file = "significant_genes_COPD_vs_Control.csv", row.names = FALSE)


# Get results and plot for Control vs COPD

plot_volcano(res_copd_df, "COPD vs Control")
head(res_copd)
rownames(dds)[1:10]  # View the first 10 row names to check for gene names

res_ra <- results(dds, contrast = c("Condition", "RA", "Control"))

# Convert the results to a data frame
res_ra_df <- as.data.frame(res_ra)

# Add gene names as a new column
res_ra_df$gene <- gene_names

# Get results and plot for Control vs RA
res_ra <- results(dds, contrast = c("Condition", "RA", "Control"))
plot_volcano(res_ra_df, "RA vs Control")
significant_genes <- res_ra_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  select(gene, log2FoldChange, padj)
print(significant_genes)
write.csv(significant_genes, file = "significant_genes_RA_vs_Control.csv", row.names = FALSE)

res_copd_ra <- results(dds, contrast = c("Condition", "RA", "Control"))

# Convert the results to a data frame
res_copd_ra_df <- as.data.frame(res_copd_ra)

# Add gene names as a new column
res_copd_ra_df$gene <- gene_names
# Get results and plot for COPD vs RA
res_copd_ra <- results(dds, contrast = c("Condition", "RA", "COPD"))
plot_volcano(res_copd_ra_df, "RA vs COPD")
significant_genes <- res_copd_ra %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  select(gene, log2FoldChange, padj)
print(significant_genes)
write.csv(significant_genes, file = "significant_genes_RA_vs_Control.csv", row.names = FALSE)

significant_genes <- res_ra_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  select(gene, log2FoldChange, padj)
print(significant_genes)
write.csv(significant_genes, file = "significant_genes_RA_vs_Control.csv", row.names = FALSE)
# Select significantly differentially expressed genes (adjust p-value < 0.05 and |log2FoldChange| > 1)
de_genes_copd <- rownames(subset(res_copd, padj < 0.05 & abs(log2FoldChange) > 1))
de_genes_ra <- rownames(subset(res_ra, padj < 0.05 & abs(log2FoldChange) > 1))

# Combine the DE gene lists for both comparisons
de_genes <- union(de_genes_copd, de_genes_ra)

# Plot heatmap of normalized counts for DE genes
heatmap.2(normalized_counts[de_genes,], 
          trace = "none", 
          Colv = TRUE, 
          scale = "row", 
          margins = c(5, 10), 
          col = bluered(75),
          main = "Heatmap of DE Genes in COPD and RA")

library(PCAtools)
p <- pca(assay(vst(dds, blind = FALSE)), metadata = colData(dds))

b1 <- biplot(p,title = "PC1 vs PC2",
             pointSize = 5 , labSize = 5 , x = "PC1" , y = "PC2")
b2 <- biplot(p,title = "PC2 vs PC3",
             pointSize = 5 , labSize = 5 , x = "PC2" , y = "PC3")
b3 <- biplot(p,title = "PC1 vs PC3",
             pointSize = 5 , labSize = 5 , x = "PC1" , y = "PC3")

b <- ggarrange(b1 , b2 , b3 , ncol = 3 , nrow = 1)



png(filename = paste0("PCAplots.png"),height = 12 , width = 25 , res = 300 , units = "in")
print(b)
dev.off()

# Step 1: List of COPD genes
copd_genes <- c(
  "SCNN1B", "CFTR", "SCNN1G", "CHRM3", "ADRB2", "SERPINA1", "CHRM1", "CHRM2", 
  "CHRM5", "CHRM4", "SCNN1A", "ADORA1", "MMP1", "MMP8", "FAM13A", "PDE3A", 
  "PDE3B", "CA4", "ADORA2A", "CA2", "CA12", "CA1", "ADORA3", "MMP7", 
  "ADRB1", "ADORA2B", "MMP13", "HMGCR", "TET2", "PDE5A", "ADGRG6", 
  "ELANE", "VDR", "TGFB2", "AGTR1", "THSD4"
)

# Step 2: List of RA genes
ra_genes <- c(
  "TNF", "IL6R", "JAK2", "TYK2", "MIF", "TRAF3IP2", "IL12B", "PTGS2", "PADI4", 
  "IL17A", "TLR7", "JAK3", "JAK1", "PTGS1", "TLR9", "CD86", "DHFR", "CD80", 
  "DHODH", "ALOX5", "MS4A1", "STAT4", "PPAT", "IL23A", "ATP4A", "FKBP1A", 
  "MC2R", "PADI2", "IL1R1", "ANKRD55", "PTPN2", "RASGRP1", "SLC6A4", "IL6", 
  "FDPS", "IRF5"
)

# Step 3: List of common genes
common_genes <- c(
  "NR3C1", "PDE4A", "PDE4D", "PDE4B", "PDE4C", "OPRM1"
)

# Step 4: Combine all gene lists
combined_genes <- unique(c(copd_genes, ra_genes, common_genes))

rownames(normalized_counts) <- gene_names

# Check the structure to confirm
str(normalized_counts)

# Print the first few row names of normalized_counts to inspect
head(rownames(normalized_counts))

# Find the common genes
common_genes <- intersect(combined_genes, rownames(normalized_counts))
print(common_genes)  # This will show you the genes that are present in normalized_counts


# Step 5: Subset the normalized count data for these genes
combined_heatmap_data <- normalized_counts[combined_genes, ] 

# Step 6: Log-transform the counts for better visualization (optional)
log_transformed_combined_data <- log2(combined_heatmap_data + 1)

sample_info <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/sample_info.csv")
# Step 7: Create sample condition annotation (for x-axis)
# Assume 'sample_info' contains a "Condition" column with RA, COPD, and Control categories
annotation_col <- data.frame(Condition = sample_info$Condition)
rownames(annotation_col) <- rownames(sample_info)  # Ensure rownames correspond to sample names

# Step 8: Define gene categories based on the combined gene list
gene_categories <- c(
  rep("COPD", length(copd_genes)),
  rep("RA", length(ra_genes)),
  rep("Common", length(common_genes))
)
names(gene_categories) <- combined_genes

# Step 9: Create a data frame for gene category annotation
annotation_row <- data.frame(Category = gene_categories[combined_genes])

# Step 10: Define custom color schemes for the annotations
gene_category_colors <- list(
  Category = c(RA = "blue", COPD = "green", Common = "purple"),  # Colors for gene categories
  Condition = c(RA = "blue", COPD = "green", Control = "grey")   # Colors for sample conditions
)

# Step 11: Define a color scale for the heatmap (e.g., from white to red)
heatmap_colors <- colorRampPalette(c("white", "red"))(50)  # 50 color shades from white to red

# Step 12: Generate the heatmap with combined genes on y-axis and sample names on x-axis
pheatmap(log_transformed_combined_data,
         cluster_rows = TRUE,    # Optionally cluster genes (y-axis)
         cluster_cols = TRUE,    # Optionally cluster samples (x-axis)
         annotation_col = annotation_col,  # Sample categories on the x-axis
         annotation_row = annotation_row,  # Gene categories on the y-axis
         annotation_colors = gene_category_colors,  # Custom colors for annotations
         color = heatmap_colors,  # Color scale from white to red
         show_rownames = TRUE,    # Show gene names on y-axis
         show_colnames = TRUE,    # Show sample names on x-axis
         fontsize_row = 8,        # Adjust font size for gene names
         fontsize_col = 10)       # Adjust font size for sample names


