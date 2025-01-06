# Load necessary libraries
library(DESeq2)
library(PCAtools)
library(ggpubr)
library(readr)

# Step 1: Read in the count data from CSV files
counts_ra <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/RA_count.csv")
counts_copd <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/COPD_case_count.csv")
counts_control <- read_csv("C:/Users/Kamalavarshini/Downloads/Deseq/Deseq/COPD_control_count.csv")

# Step 2: Prepare the data for each comparison
# Assuming the first column is "Gene" and the rest are counts
# Prepare RA vs Control
rownames(counts_ra) <- counts_ra$Gene
counts_ra <- counts_ra[, -1]  # Remove gene names

rownames(counts_control) <- counts_control$Gene
counts_control <- counts_control[, -1]

# Combine RA and Control counts
combined_ra_control <- cbind(counts_ra, counts_control)

# Create sample information for RA vs Control
sample_info_ra_control <- data.frame(
  Sample = colnames(combined_ra_control),
  Condition = c(rep("RA", ncol(counts_ra)), rep("Control", ncol(counts_control)))
)
rownames(sample_info_ra_control) <- colnames(combined_ra_control)


# Prepare COPD vs Control
rownames(counts_copd) <- counts_copd$Gene
counts_copd <- counts_copd[, -1]

# Combine COPD and Control counts
combined_copd_control <- cbind(counts_copd, counts_control)

# Create sample information for COPD vs Control
sample_info_copd_control <- data.frame(
  Sample = colnames(combined_copd_control),
  Condition = c(rep("COPD", ncol(counts_copd)), rep("Control", ncol(counts_control)))
)
rownames(sample_info_copd_control) <- colnames(combined_copd_control)
# Prepare COPD vs RA
# Combine COPD and RA counts
combined_copd_ra <- cbind(counts_copd, counts_ra)

# Create sample information for COPD vs RA
sample_info_copd_ra <- data.frame(
  Sample = colnames(combined_copd_ra),
  Condition = c(rep("COPD", ncol(counts_copd)), rep("RA", ncol(counts_ra)))
)
rownames(sample_info_copd_ra) <- colnames(combined_copd_ra)
# Step 3: Create DESeqDataSet objects for each comparison
dds_ra_control <- DESeqDataSetFromMatrix(countData = combined_ra_control,
                                         colData = sample_info_ra_control,
                                         design = ~ Condition)
head(dds_ra_control)

dds_copd_control <- DESeqDataSetFromMatrix(countData = combined_copd_control,
                                           colData = sample_info_copd_control,
                                           design = ~ Condition)

dds_copd_ra <- DESeqDataSetFromMatrix(countData = combined_copd_ra,
                                      colData = sample_info_copd_ra,
                                      design = ~ Condition)

# Step 4: Perform variance stabilization transformation for each comparison
vst_ra_control <- vst(dds_ra_control, blind = FALSE)
head(vst_ra_control)
vst_copd_control <- vst(dds_copd_control, blind = FALSE)
vst_copd_ra <- vst(dds_copd_ra, blind = FALSE)
print(colnames(vst_ra_control))
print(rownames(sample_info_ra_control))
# Step 5: Perform PCA on the transformed data for each comparison
p_ra_control <- pca(assay(vst_ra_control), metadata = sample_info_ra_control)
p_copd_control <- pca(assay(vst_copd_control), metadata = sample_info_copd_control)
p_copd_ra <- pca(assay(vst_copd_ra), metadata = sample_info_copd_ra)


# Inspect metadata
print(unique(sample_info_ra_control$Condition))
print(unique(sample_info_copd_control$Condition))
print(unique(sample_info_copd_ra$Condition))

b1 <- biplot(p_ra_control, 
             title = "RA vs Control", 
             pointSize = 5, 
             labSize = 5, 
             x = "PC1", 
             y = "PC2",
             colby = "Condition",
             )  # Color by Condition

b2 <- biplot(p_copd_control, 
             title = "COPD vs Control", 
             pointSize = 5, 
             labSize = 5, 
             x = "PC1", 
             y = "PC2",
             colby = "Condition" ,
             showLegend = TRUE)  # Color by Condition

b3 <- biplot(p_copd_ra, 
             title = "COPD vs RA", 
             pointSize = 5, 
             labSize = 5, 
             x = "PC1", 
             y = "PC2",
             colby = "Condition" , 
             showLegend = TRUE) 

# Step 7: Combine the biplots
combined_biplots <- ggarrange(b1, b2, b3, ncol = 3, nrow = 1)

# Step 8: Save combined plot to PNG
png(filename = "PCAplots.png", height = 12, width = 25, res = 300, units = "in")
print(combined_biplots)
dev.off()
