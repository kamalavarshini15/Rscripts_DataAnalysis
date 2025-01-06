BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("org.Hs.eg.db") 

install.packages("cli")
install.packages("rlang")

library(clusterProfiler)
library(org.Hs.eg.db)  # Replace with the appropriate organism package if necessary
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
#library(ggplot2) #If required


file_path <- "E:/CAG/R_scripts/genes.txt" #read.delim for Tab delimited file 
print(file_path)

output_dir <- "E:/CAG/Orobanche nana/Liver Disease"

gene_list <- readLines(file_path)

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# Convert gene symbols to Entrez IDs
gene_list <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Extract the Entrez IDs
gene_ids <- gene_list$ENTREZID

# Perform GO enrichment analysis

#Can use same code to for Molecular function and Cellular component

#orgHSegdb is used only for HOMO SAPIENS

#ont - is for the onltology we want to perform

#  ontologies <- c("BP","MF","CC")
# 
# for (ont in ontologies){
#  print(ont)
# 
# go_enrichment <- enrichGO(gene         = gene_ids,
#                           OrgDb        = org.Hs.eg.db,
#                           ont          = ont,  # Biological Process (BP), Molecular Function (MF), or Cellular Component (CC)
#                            pAdjustMethod = "BH", # Adjust p-values with Benjamini-Hochberg method
#                            pvalueCutoff = 0.05,
#                            qvalueCutoff = 0.2,
#                            readable     = TRUE)
# 
# file_name <- paste0("GO_enrichment_", ont, ".csv")
# file_path <- file.path(output_dir, file_name)
# 
# #Save results to a CSV file
# write.csv(as.data.frame(go_enrichment), file_path, row.names = FALSE)
# 
# 
# 
# # Define the file name and path
# plot_file_name <- paste0("GO_enrichment_", ont , "_barplot.jpeg")
# plot_file_path <- file.path(output_dir, plot_file_name)
# 
# 
# # Convert kegg_enrichment to a data frame (if not already)
# df <- as.data.frame(go_enrichment)
# 
# # Filter for the top 14 categories
# bar_data <- df[1:14, ]
# 
# # Create the plot with ggplot2
# p <- ggplot(bar_data, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # Flip coordinates for better readability
#   scale_fill_gradientn(colors = c("blue", "#1984c5", "#22a7f0", "#48b5c4", "#76c68f", "green")(100)) +
#   labs(title = "GO Enrichment",
#        x = ont,
#        y = "Gene Count") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Save the plot as a JPEG file
# ggsave(filename = plot_file_path, plot = p, width = 12, height = 10, dpi = 150)
# 
# # Save the plot as a JPEG file
# dp = ggsave(filename = plot_file_path, plot = p, width = 12, height = 10, dpi = 150)
# print(dp)
# 
# 
# 
# }

 
 # Define ontologies
 ontologies <- c("BP", "MF", "CC")
 
 # Loop through each ontology
 for (ont in ontologies) {
   print(ont)
   
   # Perform GO enrichment analysis
   go_enrichment <- enrichGO(gene         = gene_ids,
                             OrgDb        = org.Hs.eg.db,
                             ont          = ont,  # Biological Process (BP), Molecular Function (MF), or Cellular Component (CC)
                             pAdjustMethod = "BH", # Adjust p-values with Benjamini-Hochberg method
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2,
                             readable     = TRUE)
   
   # Save results to a CSV file
   file_name <- paste0("GO_enrichment_", ont, ".csv")
   file_path <- file.path(output_dir, file_name)
   write.csv(as.data.frame(go_enrichment), file_path, row.names = FALSE)
   
   # Define the file name and path for the bar plot
   plot_file_name <- paste0("GO_enrichment_", ont , "_barplot.jpeg")
   plot_file_path <- file.path(output_dir, plot_file_name)
   
   # Convert GO enrichment to a data frame
   df <- as.data.frame(go_enrichment)
   
   
   # Remove rows with NA values in the "Count" column (or any other relevant column)
   df_cleaned <- df[!is.na(df$Count), ]
   head(df_cleaned)
   
   
   # Filter for the top 14 categories
   bar_data <- df[1:14, ]
   
   # Create the bar plot with ggplot2 using the custom color palette
   p <- ggplot(bar_data, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
     geom_bar(stat = "identity") +
     coord_flip() +  # Flip coordinates for better readability
     scale_fill_gradientn(colors = c("#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#76c68f", "#a6d75b")) +
     labs(title = paste("GO Enrichment -", ont),
          x = ont,
          y = "Gene Count") +
     theme_minimal() + theme(axis.line = element_line(color = "black", size = 0.5),   # Line color and size for both x and y axes
                             axis.title.x = element_text(color = "black", size = 12),   # x-axis title color and size
                             axis.title.y = element_text(color = "black", size = 12),   # y-axis title color and size
                             axis.text.x = element_text(color = "black"),   # x-axis text color
                             axis.text.y = element_text(color = "black")) + 
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   # Save the plot as a JPEG file
   ggsave(filename = plot_file_path, plot = p, width = 12, height = 10, dpi = 150)
 }
 

# Perform Kegg Pathway
 


kegg_enrichment <- enrichKEGG(
  gene         = gene_ids,
  organism     = 'hsa',   # 'hsa' is the code for Homo sapiens (human)
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# View the results
head(kegg_enrichment)

# file_name <- paste0("Kegg_Pathway", ".csv")
# file_path <- file.path(output_dir, file_name)
# 
# #Save results to a CSV file
# write.csv(as.data.frame(kegg_enrichment), file_path, row.names = FALSE)
# 
# 
# 
# plot_file_name <- paste0("Kegg_Pathway","_barplot.jpeg")
# 
# plot_file_path <- file.path(output_dir, plot_file_name)
# 
# jpeg(plot_file_path, width = 1200, height = 1000, res = 150)
# 
# colors <- colorRampPalette(c("green", "yellow", "red"))(14)  # 14 is the number of categories
# 
# # Generate the bar plot with gradient colors
# bp <- barplot(kegg_enrichment, showCategory = 14, col = colors)
# 
# bp = barplot(kegg_enrichment, showCategory = 14)
# 
# #dot plot(kegg_enrichment, showCategory = 20)
# 
# #cnet plot(kegg_enrichment)
# 
# print(bp )
#
# dev.off()

# # Load required libraries
# library(clusterProfiler)
# library(ggplot2)
# library(RColorBrewer)
# 
# # Define the file name and path
# plot_file_name <- paste0("Kegg_Pathway", "_barplot.jpeg")
# plot_file_path <- file.path(output_dir, plot_file_name)
# 
# # Convert kegg_enrichment to a data frame (if not already)
# df <- as.data.frame(kegg_enrichment)
# 
# # Filter for the top 14 categories
# bar_data <- df[1:14, ]
# 
# # Create the plot with ggplot2
# p <- ggplot(bar_data, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # Flip coordinates for better readability
#   scale_fill_gradientn(colors = colorRampPalette(c("green", "blue"))(100)) +
#   labs(title = "KEGG Pathway Enrichment",
#        x = "Pathway",
#        y = "Gene Count") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Save the plot as a JPEG file
# ggsave(filename = plot_file_path, plot = p, width = 12, height = 10, dpi = 150)



# Define the file name and path
plot_file_name <- paste0("Kegg_Pathway", "_dotplot.jpeg")
plot_file_path <- file.path(output_dir, plot_file_name)

# Convert kegg_enrichment to a data frame (if not already)
df <- as.data.frame(kegg_enrichment)

# Filter for the top 20 categories
dot_data <- df[1:20, ]

# Create the dot plot with ggplot2
p <- ggplot(dot_data, aes(x = Count, y = reorder(Description, Count), size = GeneRatio, color = p.adjust)) +
  geom_point() +
  scale_fill_gradientn(colors = c("#ffb400", "#d2980d", "#a57c1b", "#786028", "#363445", "#48446e", "#5e569b", "#776bcd", "#9080ff")) +
  labs(title = "KEGG Pathway Enrichment Dot Plot",
       x = "Gene Count",
       y = "Pathway",
       color = "Adjusted p-value",
       size = "Gene Ratio") +
  theme_minimal() + theme_minimal() + theme(axis.line = element_line(color = "black", size = 0.5),   # Line color and size for both x and y axes
                                            axis.title.x = element_text(color = "black", size = 12),   # x-axis title color and size
                                            axis.title.y = element_text(color = "black", size = 12),   # y-axis title color and size
                                            axis.text.x = element_text(color = "black"),   # x-axis text color
                                            axis.text.y = element_text(color = "black")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a JPEG file
ggsave(filename = plot_file_path, plot = p, width = 12, height = 10, dpi = 150)

# Create the dot plot with ggplot2
p <- ggplot(dot_data, aes(x = Count, y = reorder(Description, Count), size = GeneRatio, color = p.adjust)) +
  geom_point() +
  scale_color_gradientn(colors = c("#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#76c68f", "#a6d75b", "#c9e52f", "#d0ee11", "#d0f400")) +
  labs(title = "KEGG Pathway Enrichment Dot Plot",
       x = "Gene Count",
       y = "Pathway",
       color = "Adjusted p-value",
       size = "Gene Ratio") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black", size = 0.5),   # Line color and size for both x and y axes
        axis.title.x = element_text(color = "black", size = 12),   # x-axis title color and size
        axis.title.y = element_text(color = "black", size = 12),   # y-axis title color and size
        axis.text.x = element_text(color = "black"),   # x-axis text color
        axis.text.y = element_text(color = "black")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a JPEG file
ggsave(filename = plot_file_path, plot = p, width = 12, height = 10, dpi = 150)


# Confirm completion
cat("GO enrichment and Kegg Pathway Analysis completed and results saved to", output_dir, "\n")









