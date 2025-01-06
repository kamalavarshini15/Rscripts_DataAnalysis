# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)

library(dplyr)   # For data manipulation
setwd("E:/CAG/mirna")
# Step 1: Load the sample names from the text file
sample_names <- readLines("input_tissue_adeno.txt")
sample_names2 <- readLines("input_plasma_adeno.txt")
#sample_names3 <- readLines("input_samples_ra.txt")
print(sample_names)

sample_names <- unlist(strsplit(sample_names, "\\s+")) 
sample_names2 <- unlist(strsplit(sample_names2, "\\s+"))
#sample_names3 <- unlist(strsplit(sample_names3, "\\s+"))
print(sample_names)

# Step 2: Load the CSV file with the sample names in the first row
# Use `header = TRUE` to read the first row as column names
data <- read.csv("counts_LAC.csv", header = TRUE)
data2 <- read.csv("counts_ADO.csv", header = TRUE)
#data3 <- read.csv("RA_count.csv", header = TRUE)
# Step 3: Extract sample names (excluding the first column)
# The first row contains the sample names starting from the second column
colnames(data)[-1] -> sample_columns
colnames(data2)[-1] -> sample_columns2
#colnames(data3)[-1] -> sample_columns3

print(sample_columns2)

print(setdiff(sample_names2, sample_columns2))
print(setdiff(sample_names3, sample_columns3))

# Step 4: Filter to keep only the relevant sample columns
# Keep the first column (identifiers) and filter for the sample names
filtered_data <- data %>%
  select(c(1, which(sample_columns %in% sample_names) + 1))
filtered_data2 <- data2 %>%
  select(c(1, which(sample_columns2 %in% sample_names2) + 1))
filtered_data3 <- data3 %>%
  select(c(1, which(sample_columns3 %in% sample_names3) + 1))

filtered_data2 <- data2 %>%
  select(c(1, which(sample_columns2 %in% sample_names2) + 1))

filtered_data3 <- data3 %>%
  select(c(1, which(sample_columns3 %in% sample_names3) + 1))

print(head(filtered_data2))
print(head(filtered_data3))

LAC_count <- filtered_data
ADO_count <- filtered_data2
#ra_count <- filtered_data3
write.csv(LAC_count,file="counts_tissue.csv")
write.csv(ADO_count,file="counts_plasma.csv")



# Merge count data (assuming first column is 'Gene' in each file)
merged_counts <- reduce(list(LAC_count, ADO_count), full_join, by = "miRNA_ID")
# set gene column to rownames
merged_counts <- merged_counts %>%column_to_rownames(var = "miRNA_ID")
colnames(merged_counts) <- gsub("_gene_counts", "", colnames(merged_counts))
write.csv(merged_counts,file="Deseq_dataset_tissue_vs_plasma.csv")
# Load sample information
sample_info <- read_csv("sample_info_NEW.csv")  # Load your sample information CSV
combined_sample_names <- c( sample_names , sample_names2)
print(combined_sample_names)

sample_info_filtered <- sample_info %>% filter(Sample %in% combined_sample_names)
print(sample_info_filtered)

rownames(sample_info_filtered) <- sample_info_filtered$Sample

# Ensure columns in counts match samples in sample_info
counts <- merged_counts[,rownames(sample_info_filtered)]

# Create DESeq2 dataset using sample condition from sample_info
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info_filtered, design = ~ Condition)

# Run the DESeq2 pipeline
dds <- DESeq(dds)
head(dds)

# Extract normalized counts (optional step, DESeq2 internally uses normalized data for analysis)
normalized_counts <- counts(dds, normalized = TRUE)

# Print normalized counts (optional, just to see how they look)
head(normalized_counts)

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

res_tissue_plasma <- results(dds, contrast = c("Condition", "ADO_tissue", "ADO_plasma"))

res_tissue_plasma_df <- as.data.frame(res_tissue_plasma)
mirna_names <- merged_counts$miRNA_ID

# Add gene names as a new column
res_tissue_plasma_df$miRNA_ID <- mirna_names

plot_volcano(res_tissue_plasma_df, "Adenocarcinoma:Tissue vs Plasma")



  # Extract rownames, which are the gene names



# Plotting function
plot_volcano <- function(res, title) {
  res$logp <- -log10(res$pvalue)
  res$significance <- ifelse(res$padj < 0.05 & res$log2FoldChange > 2, "Upregulated",
                             ifelse(res$padj < 0.05 & res$log2FoldChange < -2, "Downregulated", "Not Significant"))
  
  # Create a volcano plot with color coding for upregulated and downregulated genes
  p <- ggplot(res, aes(x = log2FoldChange, y = logp, color = significance)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal()
  # Extract significant genes based on log2FoldChange and padj thresholds
 
  
  # Add gene names to significant points (optional)
  # Uncomment if you want to label significant genes
  # p + geom_text_repel(data = subset(res, significance != "Not Significant"),
  #                     aes(label = gene), 
  #                     size = 3, 
  #                     max.overlaps = Inf)  # Avoid overlapping labels
  
  return(p)
}

# Generate the volcano plot
volcano_plot <- plot_volcano(res_tissue_plasma_df, "NLAC vs LAC")
print(volcano_plot)

filtered_genes <- subset(res_tissue_plasma_df, padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2))
write.csv(filtered_genes, file="Significant mirna_adeno_tissue_plasma.csv")

# Add a column for upregulated or downregulated based on log2FoldChange
filtered_genes$regulation_status <- ifelse(filtered_genes$log2FoldChange > 0, "Upregulated", 
                                           ifelse(filtered_genes$log2FoldChange < 0, "Downregulated", "No Change"))

write.csv(filtered_genes, file="Significant mirna_adeno_tissue_plasma_regulation.csv")
res_tissue_plasma_df$condition 
<- ifelse(grepl("Tissue", rownames(res_tissue_plasma_df)), "ADO_tissue", "ADO_plasma")

#

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


# Step 1: Read the data from a CSV file
# Replace 'your_file.csv' with the actual file path
mirnet_data <- read.csv("info.csv")

# Step 2: Apply the comparison logic
# Assuming your CSV file has columns 'mean_disease_A' and 'mean_disease_B'
mirnet_data$status <- ifelse(mirnet_data$Tissue > mirnet_data$Plasma, 
                      "upregulated for Tissue downregulated for Plasma", 
                      "upregulated for Plasma downregulated for Tissue")

# Step 3: View the updated data
print(mirnet_data)

# Optional Step 4: Write the updated data back to a CSV file
# Replace 'updated_file.csv' with your desired output file path
write.csv(mirnet_data, "updated_file.csv", row.names = FALSE)
# Read the files
file1 <- read.csv("Significant mirna_adeno_tissue_plasma_regulation.csv")
file2 <- read.csv("updated_file.csv")
# Select the necessary columns: miRNA_ID and log2 fold
file1_selected <- file1[, c("miRNA_ID", "log2FoldChange")]
# Merge the two data frames by miRNA_ID
merged_data <- merge(file2, file1_selected, by = "miRNA_ID", all.x = TRUE)
# Write the merged data to a new CSV file
write.csv(merged_data, "merged_file.csv", row.names = FALSE)
