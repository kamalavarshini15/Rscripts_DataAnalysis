
setwd("E:/CAG/Variant_Disease_Prediction_AI")



# Load libraries
library(readxl)
library(writexl)
library(dplyr)

# Read the sheets into data frames
sheet1 <- read_excel("Parameters.xlsx", sheet = "Sheet1")
sheet2 <- read_excel("Parameters.xlsx", sheet = "Sheet2")

# Add row numbers to Sheet1
sheet1 <- sheet1 %>%
  mutate(RowNumber_Sheet1 = row_number())

# Filter rows where the data starts with ""
filtered_data <- sheet1 %>%
  filter(grepl("^GERPNR=", GERPNR)) 


# Map filtered data to Sheet2 based on the 'ID' column
mapped_data <- sheet2 %>%
  left_join(filtered_data, by = "Num")

# Write the updated data to an Excel file
write_xlsx(mapped_data, "updated_file.xlsx")

# View the first few rows of the mapped data
head(mapped_data)

