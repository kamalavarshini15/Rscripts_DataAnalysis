# Install necessary packages (if not already installed)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

# Load the packages
library(readxl)
library(dplyr)

# Define file paths for the two Excel files
file1 <- "HC.xlsx"  # Replace with the path to your first Excel file
file2 <- "variant_info_new.xlsx"  # Replace with the path to your second Excel file

# Read the sheets from the Excel files
data1 <- readxl::read_excel(file1)  # Replace '1' with sheet name or index if needed
data2 <- readxl::read_excel(file2)

# Combine the two data frames based on IDs
# Assuming the column with IDs is named 'ID' in both files
combined_data <- data1 %>%
  inner_join(data2, by = "dbSNP ID")  # Use other joins like left_join, right_join, or full_join if needed

# View the combined data
print(combined_data)

# Save the combined data to a new Excel file
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
library(writexl)

write_xlsx(combined_data, "combined_data.xlsx")  # Replace with your desired output file name
