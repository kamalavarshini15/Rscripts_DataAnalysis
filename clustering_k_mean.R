# Load required libraries
library(vcfR)
library(tidyverse)
library(cluster)
library(factoextra)
setwd("E:/CAG/FRGA_Dataset_first_10/GORA_DHAN_2IRGC_66269-1")
# Step 1: Load VCF files
vcf_files <- list.files(path = "E:/CAG/FRGA_Dataset_first_10/GORA_DHAN_2IRGC_66269-1", pattern = "_SNPs.*\\.vcf$", full.names = TRUE)
print(vcf_files)
# Step 2: Function to convert genotype data into numeric format
convert_genotype <- function(gt) {
  if (gt == "0/0") {
    return(0)  # Homozygous reference
  } else if (gt == "1/1" || gt == "2/2" || gt == "1|1" || gt == "2|2") {
    return(2)  # Homozygous variant
  } else if (gt == "0/1" || gt == "1/0" || gt == "1/2" || gt == "1|2" || gt == "0|1" || gt == "1|0") {
    return(1)  # Heterozygous variant
  } else {
    return(NA)  # Missing data
  }
}

# Step 3: Extract genotype data from all VCF files and convert to numeric format
genotypes_list <- lapply(vcf_files, function(file) {
  vcf_data <- read.vcfR(file)
  genotypes <- extract.gt(vcf_data)
  genotypes_numeric <- apply(genotypes, c(1, 2), function(x) convert_genotype(x))
  return(genotypes_numeric)
})

max_rows <- max(sapply(genotypes_list, nrow))


#fill missing places with NA 
genotypes_list <- lapply(genotypes_list, function(x) {
  if (nrow(x) < max_rows) {
    x <- rbind(x, matrix(NA, nrow = max_rows - nrow(x), ncol = ncol(x)))
  }
  return(x)
})

genotypes_combined <- do.call(cbind, genotypes_list)



# Step 5: Scale the genotype matrix (important for K-means)
genotypes_scaled <- scale(genotypes_combined)

# Step 6: Perform K-means clustering 
set.seed(123)  # For reproducibility

genotypes_scaled[is.na(genotypes_scaled)] <- apply(genotypes_scaled, 2, function(x) mean(x, na.rm = TRUE))

kmeans_result <- kmeans(genotypes_scaled, centers = 3)  # You can change the number of clusters (centers) as needed

# Step 7: View the clustering result (which sample belongs to which cluster)
kmeans_result$cluster  # Cluster assignments for each SNP

# Step 8: Combine the genotype data with the cluster information
genotypes_with_clusters <- data.frame(genotypes_scaled, cluster = kmeans_result$cluster)

# Step 9: View the first few rows of the data with cluster assignments
head(genotypes_with_clusters)

# Assuming kmeans_result is already created
sil_score <- silhouette(kmeans_result$cluster, dist(genotypes_scaled))
fviz_silhouette(sil_score)
avg_sil_score <- mean(sil_score[, 3])  # The third column contains the silhouette width
print(paste("Average Silhouette Score: ", avg_sil_score))




# Assuming pca_data is already created
write.csv(genotypes_with_clusters, "genotypes_with_clusters_GORA_DHAN_2IRGC_66269-1.csv", row.names = TRUE)


# Step 10: Optional - Plot the results using PCA for dimensionality reduction (2D plot)
pca_result <- prcomp(genotypes_scaled)  # Perform PCA
pca_data <- data.frame(pca_result$x, cluster = factor(kmeans_result$cluster))


# Plot PCA with cluster coloring
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "PCA of SNP Data with K-means Clustering")

# Step 11: Optional - Determine the optimal number of clusters using the elbow method
wss <- sapply(1:10, function(k) {
  kmeans(genotypes_scaled, centers = k)$tot.withinss
})

# Plot the elbow curve to determine the optimal number of clusters
plot(1:10, wss, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters", ylab = "Total within-cluster sum of squares")
