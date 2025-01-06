# Load necessary libraries
library(DESeq2)
library(PCAtools)
library(ggpubr)
library(readr)
setwd("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq")
# Step 1: Read in the count data from CSV files
counts_ra <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/RA_count.csv")
counts_copd <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/COPD_case_count.csv")
counts_control <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/COPD_control_count.csv")

# Step 2: Prepare the data for each comparison
# Assuming the first column is "Gene" and the rest are counts
# Set gene names as row names and remove the first column (gene names)
rownames(counts_ra) <- counts_ra$Gene
counts_ra <- counts_ra[, -1]

rownames(counts_copd) <- counts_copd$Gene
counts_copd <- counts_copd[, -1]

rownames(counts_control) <- counts_control$Gene
counts_control <- counts_control[, -1]

# Step 3: Combine COPD, RA, and Control counts into a single matrix
combined_copd_ra_control <- cbind(counts_copd, counts_ra, counts_control)

# Step 4: Create sample information for COPD, RA, and Control
sample_info_copd_ra_control <- data.frame(
  Sample = colnames(combined_copd_ra_control),
  Condition = c(rep("COPD", ncol(counts_copd)), 
                rep("RA", ncol(counts_ra)), 
                rep("Control", ncol(counts_control)))
)
rownames(sample_info_copd_ra_control) <- colnames(combined_copd_ra_control)

# Step 5: Create DESeqDataSet object for all three conditions
dds_copd_ra_control <- DESeqDataSetFromMatrix(countData = combined_copd_ra_control,
                                              colData = sample_info_copd_ra_control,
                                              design = ~ Condition)

# Step 6: Perform variance stabilizing transformation (VST)
vst_copd_ra_control <- vst(dds_copd_ra_control, blind = FALSE)

# Step 7: Perform PCA on the transformed data
p_copd_ra_control <- pca(assay(vst_copd_ra_control), metadata = sample_info_copd_ra_control)

# Step 8: Create a PCA biplot colored by condition (COPD, RA, and Control)
b_copd_ra_control <- biplot(p_copd_ra_control, 
                            title = "COPD vs RA vs Control", 
                            pointSize = 5, 
                            labSize = 5, 
                            x = "PC1", 
                            y = "PC2",
                            colby = "Condition",  # Color by Condition
                            legendPosition = "right") 
                            

# Step 9: Save the PCA plot as a PNG file
png(filename = "PCAplot_COPD_RA_Control.png", height = 12, width = 15, res = 300, units = "in")
print(b_copd_ra_control)
dev.off()
