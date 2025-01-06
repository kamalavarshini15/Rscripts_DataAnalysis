# Load necessary libraries
library(rentrez)
library(XML)

# 1. Search for Taxonomy Record of 'Tetrahymena thermophila'
# We are looking for the species rank for 'Tetrahymena thermophila'
Tt <- entrez_search(db = "taxonomy", term = "(Methanothrix soehngenii[ORGN]) AND Species[RANK]")

# Check if search returned results
if (Tt$count > 0) {
  
  # 2. Fetch the Taxonomy Record in XML format
  tax_rec <- entrez_fetch(db = "taxonomy", id = Tt$ids, rettype = "xml", parsed = TRUE)
  
  # Print class of tax_rec to verify it's parsed XML data
  print(class(tax_rec))
  
  # 3. Convert the XML result to a list for easier manipulation
  tax_list <- xmlToList(tax_rec)
  
  # Extract the genetic code for the species if available
  genetic_code <- tax_list$Taxon$GeneticCode
  
  # Handle the genetic code printing appropriately
  if (!is.null(genetic_code)) {
    cat("Genetic Code:\n")
    print(genetic_code)  # Use print to handle lists
  } else {
    cat("No Genetic Code information found.\n")
  }
  
  # 4. Extract lineage information using XPath
  # Get the scientific names from the lineage
  tt_lineage <- xpathSApply(tax_rec, "//LineageEx/Taxon/ScientificName", xmlValue)
  
  # Display the first few entries of the lineage
  cat("Lineage (first 4 entries):\n")
  print(tt_lineage[1:4])  # Show first 4 entries of the lineage
  
  # 5. Print the full lineage (optional, if you want to see the entire path)
  cat("\nFull Lineage:\n")
  print(tt_lineage)
  
} else {
  # No records found
  cat("No taxonomy records found for 'Tetrahymena thermophila'\n")
}
print(tt_lineage)
