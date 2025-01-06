install.packages("rentrez")


# Load necessary libraries
library(rentrez)

# Load DSMZ data from an Excel file (adjust the path as necessary)
dsmz_data <- read.csv("E:/CAG/16s_Database/Full_Taxonomy_Data.csv")  

# Create a directory to store the downloaded 16S rRNA sequences
output_dir <- "E:/CAG/16s_Database"
dir.create(output_dir, showWarnings = FALSE)

# Function to download 16S rRNA sequences for a given genus and species
download_16S_rRNA <- function(genus, species = NULL) {
  # Construct the search term
  search_query <- genus
  if (!is.null(species) && species != "") {
    search_query <- paste(genus, species, "16S","complete sequence","ribosomal RNA")
  }
  
  message(paste("Searching for:", search_query))
  
  # Search for the sequences in the nucleotide database
  search_results <- tryCatch({
    entrez_search(db = "nucleotide", term = search_query, retmax = 10)
  }, error = function(e) {
    message("Error during search: ", e$message)
    return(NULL)
  })
  
  # Debugging output to inspect the search results structure
  if (!is.null(search_results)) {
    print(search_results)
  }
  
  # Check if the search results are valid
  if (!is.null(search_results) && !is.null(search_results$ids) && length(search_results$ids) > 0) {
    # Fetch the sequences using the IDs from the search results
    sequences <- tryCatch({
      entrez_fetch(db = "nucleotide", id = search_results$ids, rettype = "fasta", retmode = "text")
    }, error = function(e) {
      message("Error during fetch: ", e$message)
      return(NULL)
    })
    
    # Save the sequences to a FASTA file
    if (!is.null(sequences)) {
      fasta_file <- file.path(output_dir, paste0(genus, "_", ifelse(!is.null(species) && species != "", species, "unknown"), "_16S.fasta"))
      writeLines(sequences, fasta_file)
      message(paste("Downloaded 16S rRNA sequences for:", genus, species))
    } else {
      message(paste("No sequences returned for:", genus, species))
    }
  } else {
    message(paste("No sequences found for:", genus, species))
  }
}



# Iterate over rows in the DSMZ dataset and extract genus and species
for (i in 1:nrow(dsmz_data)) {
  genus <- dsmz_data$genus_name[i]
  species <- dsmz_data$sp_epithet[i]
  
  # Check if genus is not NA
  if (!is.na(genus) && genus != "") {
    download_16S_rRNA(genus, species)
  } else {
    message("No valid genus name found at row: ", i)
  }
}

