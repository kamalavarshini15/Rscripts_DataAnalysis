# Function to combine all FASTA files in a directory into a single file
combine_fasta_files <- function(input_directory, output_file) {
  # Get all FASTA files from the directory
  fasta_files <- list.files(path = input_directory, pattern = "\\.fasta$", full.names = TRUE)
  
  # Check if any FASTA files were found
  if (length(fasta_files) == 0) {
    message("No FASTA files found in the directory.")
    return()
  }
  
  # Open the output file for writing
  output_connection <- file(output_file, open = "w")
  
  # Loop through each FASTA file and append to the output file
  for (fasta_file in fasta_files) {
    message(paste("Processing file:", fasta_file))
    
    # Read the FASTA file content
    fasta_content <- readLines(fasta_file)
    
    # Write the content to the output file
    writeLines(fasta_content, output_connection)
    
    # Add a newline between the contents of different FASTA files
    writeLines("", output_connection)
  }
  
  # Close the output file
  close(output_connection)
  
  message(paste("All FASTA files combined into:", output_file))
}

# Example usage
input_directory <- "E:/CAG/16s_Database/test_data"  # Specify the directory containing FASTA files
output_file <- "combined_sequences.fasta"      # Specify the output file name
combine_fasta_files(input_directory, output_file)
