install.packages("recipes")  # Dependencies for "Caret" package 
install.packages("rQSAR")
install.packages("dplyr")
install.packages("corrplot")
install.packages("tibble")
install.packages("gridExtra")
install.packages("rcdk")
install.packages("ggplot2")
install.packages("caret")
install.packages("pls")
install.packages("randomForest")
install.packages("car")



library(dplyr)
library(corrplot)
library(tibble)
library(gridExtra)
library(rJava)
library(rcdklibs)
library(ChemmineR) #(Only if need)
library(rcdk)
library(ggplot2)
library(caret)  #(Only if need)
library(pls)
library(randomForest)
library(rQSAR)
library(FSelector)


output_dir <- "E:/CAG/R_scripts/qsar" #Output Directory 

#DATA PREPARATION 

#(Canonical smiles were parsed into molecule objects and calculated chemical descriptors using "rcdk" )
# Load CSV file containing molecule names and canonical SMILES
input_file <- "E:/CAG/mol1.csv"  # Replace with the path to your CSV file
molecule_data <- read.csv(input_file)
print(molecule_data)

#If the Input are given in a text box as a string 

# smiles_list <- c(
#   "CCO",    # Ethanol
#   "CC(C)O", # Isopropanol
#   "CCC(=O)O", # Propionic acid
#   "CCN",    # Ethylamine
#   "C1CC1",  # Cyclopropane
#   "C1=CC=CC=C1", # Benzene
#   "CC(=O)O", # Acetic acid
#   "C=O",    # Formaldehyde
#   "CC#N",   # Acetonitrile
#   "O=C=O"   # Carbon dioxide
# )

#molecules <- parse.smiles(smiles_list)

# Parse SMILES to molecule objects
molecules <- parse.smiles(molecule_data$Canonical.SMILES)
# Filter out any NULL entries 
molecules <- Filter(Negate(is.null), molecules)
# Check if molecules are correctly parsed
if (length(molecules) == 0) {
  stop("No valid molecules were parsed from the SMILES strings.")
}
# Print the valid molecules
print(molecules)

# List available descriptors
all_desc_names <- get.desc.names()
# print(all_desc_names)  # View all available descriptors

# # Check specific descriptor names
# required_desc_names <- c(
#   "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",              
#   "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor",
#   "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor"
#   
# )

# # Ensure all required descriptors are available
# valid_desc_names <- intersect(required_desc_names, all_desc_names)
# print(valid_desc_names)
# # Calculate descriptors for valid molecules
# descriptors <- eval.desc(molecules, valid_desc_names)

# Calculate all descriptors for the valid  molecules
descriptors <- eval.desc(molecules, all_desc_names)
#Save descriptors as CSV file
file_name <- paste0("discriptors", ".csv")
file_path <- file.path(output_dir, file_name)
write.csv(as.data.frame(descriptors), file_path, row.names = FALSE)
# Remove columns with all NA values
descriptors <- descriptors[, colSums(!is.na(descriptors)) > 0]
# Ensure all columns are numeric
descriptors <- descriptors[, sapply(descriptors, is.numeric)]
# Filter out columns with all zero values using dplyr
filtered_descriptors <- descriptors %>%
  select(where(~ any(. != 0)))
# Check the filtered data
print(filtered_descriptors)
write.csv(filtered_descriptors, "data_filtered.csv", row.names = TRUE)

# Example activity data (Random Values used)
# activity_data <- c(1.2, 3.4, 2.1, 4.5, 3.3, 8.7, 9.6, 5.45, 6.24 , 8.3)  # Replace with your actual activity data
activity_data <- (molecule_data$Activity)
print(activity_data)

#DATA PREPROCESSING 
#( Converted Descriptors into Data frame and combined activity data with structural data as CSV for Model building )

# Prepare descriptor data frame
desc_df_list <- lapply(filtered_descriptors, function(desc) {
  # Ensure that each descriptor has the same length
  as.data.frame(t(desc), check.names = TRUE)
})
# Combine all descriptor data frames
desc_df <- do.call(rbind, desc_df_list)
# Convert activity_data to a data frame row
activity_row <- as.data.frame(t(activity_data))
# Add row names for clarity, if needed
rownames(activity_row) <- "Activity"
# Append activity_row to desc_df
data <- rbind(desc_df, activity_row)
# Print the new data frame
print(data)
#Took transpose
t_data <- t(data)
print(t_data)

t_data_df <- as.data.frame(t_data)

# Check the structure of the data frame
str(t_data)

# Check the column names
print(names(t_data))

# Verify that the "Activity" column exists
if ("Activity" %in% names(t_data)) {
  # Exclude the "Activity" column and get the attributes
  attributes <- names(t_data)[-which(names(t_data) == "Activity")]
  print(attributes)
} else {
  print("Activity column not found in the data frame.")
}





# Check if attributes are specified correctly
attributes <- names(t_data_df)[-which(names(t_data_df) == "Activity")]

print(attributes)

# Define the evaluation function (e.g., correlation with Activity)
eval.fun <- function(attributes) {
  selected_data <- t_data_df[, c("Activity", attributes)]
  cor(selected_data$Activity, selected_data[, 2])
}

print(eval.fun)




attributes <- names(t_data_df)[-which(names(t_data_df) == "Activity")]
print(attributes) 


best_descriptors <- best.first.search(
  attributes = attributes,
  eval.fun = eval.fun,
  max.backtracks = 5
)


eval.fun <- function(attrs) {
  # Fit a linear model using the selected attributes
  formula <- as.formula(paste("Activity ~", paste(attrs, collapse = " + ")))
  model <- lm(formula, data = t_data_df)
  
  # Return the R-squared value as the evaluation metric
  return(summary(model)$r.squared)
}



# Apply best.first.search to select descriptors
best_descriptors <- best.first.search(
  attributes = names(t_data_df)[-1],  # Exclude the target variable "Activity"
  eval.fun = eval.fun,
  max.backtracks = 5  # Limit the number of backtracks
)

# Print the selected descriptors
print(best_descriptors)


# # Compute the correlation matrix
# cor_matrix <- cor(t_data, use = "complete.obs")
# is.matrix(cor_matrix)  # Should return TRUE
# # Example: Remove the row and column corresponding to "MDEC.11"
# cor_matrix <- cor_matrix[!(rownames(cor_matrix) %in% "MDEC.11"), 
#                          !(colnames(cor_matrix) %in% "MDEC.11")]
# 
# cor_matrix <- cor_matrix[!(rownames(cor_matrix) %in% "khs.sssCH"), 
#                          !(colnames(cor_matrix) %in% "khs.sssCH")]
# 
# # #Want to check relation you want 
# # activity_correlation <- cor_matrix["Activity","ALogP" ]
# # # Print the correlation matrix
# # print(activity_correlation)
# 
# # Strep 1: t_data is your data frame containing the descriptors
# cor_matrix <- cor(t_data)
# # Step 2: Identify highly correlated pairs (e.g., |R| > 0.8)
# high_corr_pairs <- which(abs(cor_matrix) > 0.8 , arr.ind = TRUE) 
# # Remove diagonal elements to avoid self-correlation
# high_corr_pairs <- high_corr_pairs[high_corr_pairs[,1] != high_corr_pairs[,2], ]
# 
# #Perform STEP 3 and 4 for High Correlation or else do STEP 5 and 6 for Low Correlation
# #PCA 
# 
# # Step 3: Extract unique descriptors involved in high correlations
# # Use a set to keep track of descriptors to include
# high_corr_descriptors <- unique(c(high_corr_pairs[,1], high_corr_pairs[,2]))
# # Step 4: Create a subset of descriptors with high inter-descriptor correlations
# high_corr_descs <- t_data[, high_corr_descriptors]
# 
# 
# # Step 5: Identify descriptors to remove (keep only one from each highly correlated pair)
# descriptors_to_remove <- unique(high_corr_pairs[,2])
# # Step 6: Create a subset of descriptors with low inter-descriptor correlations
# low_corr_descs <- t_data[, -descriptors_to_remove]
# 
# 
# # Generate the heat map
# heatmap(low_corr_descs, 
#         main = "Subset Correlation Matrix Heatmap", 
#         col = colorRampPalette(c("blue", "white", "red"))(101), 
#         scale = "none", 
#         margins = c(10, 10))
# print(heatmap)




set.seed(123)  # For reproducibility
train_index <- createDataPartition(t_data_df$Activity, p = 0.8, list = FALSE)
train_data <- t_data_df[train_index, ]
test_data <- t_data_df[-train_index, ]



# Assuming best_descriptors contains the selected descriptors
# Make sure the list of descriptors doesn't include "Activity"
best_descriptors <- best_descriptors[best_descriptors != "Activity"]


# Example: Identifying and removing highly collinear features
corr_matrix <- cor(train_data_best[, -1])
high_corr <- findCorrelation(corr_matrix, cutoff = 0.9)
train_data_best <- train_data_best[, -high_corr]


# Subset train_data_cleaned and test_data to include only the selected descriptors plus "Activity"
train_data_best <- train_data[, c("Activity", best_descriptors)]
print(train_data_best)
test_data_best <- test_data[, c("Activity", best_descriptors)]
print(test_data_best)

# Build the MLR model using the selected descriptors
lm_model_best <- lm(Activity ~ ., data = train_data_best)

# Print the summary of the model to check the coefficients and statistics
summary(lm_model_best)

# Predict on the test set using the model
predictions <- predict(lm_model_best, newdata = test_data_best)

# Print the predictions to check the model's performance
print(predictions)

# Optional: Evaluate model performance with R-squared or other metrics
r_squared <- summary(lm_model_best)$r.squared
cat("R-squared of the model on training data:", r_squared, "\n")

# Calculate R-squared for test data
test_r_squared <- cor(test_data_best$Activity, predictions)^2
cat("R-squared on test data:", test_r_squared, "\n")


# Step 1: Extract the intercept and coefficients
intercept <- coef(lm_model_best)[1]
coefficients <- coef(lm_model_best)[-1]  # All coefficients except the intercept

# Step 2: Calculate the predicted values manually
manual_predictions <- intercept + as.matrix(test_data_best[, -1]) %*% coefficients

# Step 3: Compare the manually calculated predictions with model predictions
model_predictions <- predict(lm_model_best, newdata = test_data_best)

# Print both to check if they match
print(manual_predictions)
print(model_predictions)

# Step 4: Validate model fit by calculating errors
error <- model_predictions - test_data_best$Activity

# Calculate Mean Squared Error (MSE)
mse <- mean(error^2)
cat("Mean Squared Error:", mse, "\n")

# Optional: Visualize the actual vs. predicted values
plot(test_data_best$Activity, model_predictions, main = "Actual vs Predicted Values",
     xlab = "Actual Activity", ylab = "Predicted Activity", pch = 19)
abline(0, 1, col = "red")  # Add a y=x line for reference





















# Load necessary library
library(caret)  # For data splitting


# low_corr_descs <- as.data.frame(low_corr_descs)
# data_model <- low_corr_descs

low_corr_descs <- as.data.frame(low_corr_descs)
data_model <- low_corr_descs



matching_rows <- rownames(low_corr_descs)  # Example, if rownames are used
# Extract 'Activity' for the corresponding rows
activity_column <- t_data[matching_rows, "Activity"]
# Add 'Activity' column to 'low_corr_descs'
data_model <- cbind(Activity = activity_column, low_corr_descs)
# Step 6: Split the data into training and testing sets (80% train, 20% test)
set.seed(123)  # For reproducibility
train_index <- createDataPartition(data_model$Activity, p = 0.8, list = FALSE)
train_data <- data_model[train_index, ]
test_data <- data_model[-train_index, ]


dim(train_data)
dim(test_data)

# Example: Removing aliased variables manually
train_data_cleaned <- train_data[, !names(train_data) %in% c("khs.dNH", "khs.ssS")]
print(train_data_cleaned)
# Remove a column by name
names(train_data_cleaned)[names(train_data_cleaned) == "Activity.1"] <- "Activity"
print(train_data_cleaned)

lm_model <- lm(Activity ~ ., data = train_data_cleaned)

print(lm_model)

# Step 7: Build the MLR model using the training data
# mlr_model <- lm(Activity ~ ., data = train_data)


# Step 8: Predict on the test data
predictions <- predict(lm_model, newdata = test_data)

# Step 9: Evaluate the model
# Calculate R-squared and RMSE for the test set
r_squared <- cor(test_data$Activity, predictions)^2
rmse <- sqrt(mean((test_data$Activity - predictions)^2))

# Print the summary of the model and evaluation metrics
summary(mlr_model)
cat("Test R-squared:", r_squared, "\n")
cat("Test RMSE:", rmse, "\n")






# Fit the PLS model on the training data
pls_model <- plsr(Activity ~ ., data = train_data, ncomp = 2)  # Adjust ncomp as needed

# Make predictions on the test data
predictions <- predict(pls_model, newdata = test_data, ncomp = 2)

# Evaluate the model
actuals <- test_data$Activity
mse <- mean((actuals - predictions)^2)
rmse <- sqrt(mse)

# Print evaluation metrics
cat("RMSE:", rmse, "\n")






# Fit the Random Forest model on the training data

train_data_clean <- na.omit(train_data)
rf_model <- randomForest(Activity ~ ., data = train_data_clean, ntree = 100)


# Make predictions on the test data
predictions <- predict(rf_model, newdata = test_data)

# Evaluate the model
actuals <- test_data$Activity
mse <- mean((actuals - predictions)^2)
rmse <- sqrt(mse)
r_squared <- cor(actuals, predictions)^2

# Print evaluation metrics
cat("RMSE:", rmse, "\n")
cat("R-squared:", r_squared, "\n")



































 






















































