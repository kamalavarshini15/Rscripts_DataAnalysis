# Rscripts_DataAnalysis
This repo containts the r scripts for Data analysis for Genomic , proteomics problems  etc 

# R Scripts Repository
This repository contains a collection of R scripts developed for various bioinformatics, statistical analysis, and data visualization tasks. Each script is tailored for specific analyses and workflows, providing robust solutions for biological and chemical data processing.

# List of Scripts and Descriptions

**1. 16s_dmsz_ncbi.R**
- Purpose: Retrieves and processes 16S rRNA sequence data from DSMZ and NCBI databases.
- Features:
Fetches 16S rRNA sequences by ID.
Formats sequences for downstream analysis.
- Dependencies: httr, jsonlite, seqinr.
- Use Case: Constructing reference databases or analyzing microbial communities.

**2. Deseq_COPD_RA_RNAseq.R**
Purpose: Performs differential expression analysis on RNA-seq data from COPD and RA studies.
Features:
Uses DESeq2 for identifying differentially expressed genes.
Includes data preprocessing and normalization steps.
Dependencies: DESeq2, ggplot2, pheatmap.
Use Case: Gene expression analysis for COPD and RA studies.

**3. GO_Analysis_Rcode.R**
Purpose: Conducts Gene Ontology (GO) enrichment analysis.
Features:
Analyzes functional categories of genes.
Visualizes GO terms using bar plots and dot plots.
Dependencies: clusterProfiler, org.Hs.eg.db, ggplot2.
Use Case: Functional annotation of gene sets.

**4. PCAplot_combined.R**
Purpose: Generates combined PCA plots for visualizing data clusters.
Features:
Combines multiple datasets for PCA visualization.
Customizable color schemes and labels.
Dependencies: ggplot2, FactoMineR, factoextra.
Use Case: Exploratory data analysis for multi-condition experiments.

**5. QSAR_Model_Copy.R**
Purpose: Builds and validates QSAR models.
Features:
Implements MLR, PLS, Ridge, and Lasso regression.
Handles descriptor data preprocessing.
Dependencies: glmnet, caret, rcdk.
Use Case: Chemical descriptor modeling and activity prediction.

**6. clustering_k_mean.R**
Purpose: Performs k-means clustering on biological or chemical datasets.
Features:
Optimal number of clusters using the elbow method.
Visualizes clusters using scatter plots.
Dependencies: stats, factoextra, ggplot2.
Use Case: Grouping data points based on similarity.

**7. heatmap_script.R**
Purpose: Creates heatmaps for gene expression or descriptor data.
Features:
Color-coded representation of data matrices.
Customizable clustering options.
Dependencies: pheatmap, ComplexHeatmap, ggplot2.
Use Case: Visualizing patterns in expression data.

**8. indel_circos_plot.R**
Purpose: Generates circos plots for insertion-deletion (INDEL) data.
Features:
Visualizes genomic variation data in circular layout.
Supports customization of tracks and labels.
Dependencies: circlize, GenomicRanges.
Use Case: Genomic variation visualization.

**9. lpsn_retrival.R**
Purpose: Retrieves taxonomic data from the List of Prokaryotic Names with Standing in Nomenclature (LPSN).
Features:
Fetches and formats prokaryotic taxonomic information.
Dependencies: httr, jsonlite.
Use Case: Taxonomic data integration.

**10. pca_plot_script.R**
Purpose: Generates PCA plots for single datasets.
Features:
Visualizes principal components with variance explained.
Dependencies: ggplot2, FactoMineR, factoextra.
Use Case: Dimensionality reduction and visualization.

**11. r_Script_combinetwoexcel_by_IDs.R**
Purpose: Merges two Excel files by common IDs.
Features:
Matches and combines data from two spreadsheets.
Outputs a combined Excel file.
Dependencies: readxl, writexl, dplyr.
Use Case: Data integration and cleanup.

**12. r_script_model_new.R**
Purpose: Builds predictive models using advanced statistical methods.
Features:
Implements and compares models like Random Forest, SVM, and MLR.
Dependencies: caret, randomForest, e1071.
Use Case: Predictive modeling for biological or chemical datasets.

**13. retrive_fasta_sequences.R**
Purpose: Retrieves FASTA sequences from online databases or local files.
Features:
Fetches sequences by ID.
Formats and saves sequences for analysis.
Dependencies: httr, seqinr.
Use Case: Sequence data retrieval.

**14. script_r_filtering_excel.R**
Purpose: Filters Excel data based on specified criteria.
Features:
Applies custom filters to large Excel files.
Dependencies: readxl, dplyr, writexl.
Use Case: Preprocessing Excel datasets.

**15. snp_circos_plot.R**
Purpose: Creates circos plots for SNP data.
Features:
Visualizes SNPs in a circular genomic context.
Dependencies: circlize, GenomicRanges.
Use Case: SNP data visualization.

**16. taxonomy_retrival.R**
Purpose: Automates taxonomy retrieval for biological sequences.
Features:
Fetches taxonomic details from online or local sources.
Dependencies: httr, jsonlite, dplyr.
Use Case: Taxonomic annotation of sequence data.

