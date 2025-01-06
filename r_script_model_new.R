# Load required libraries
library(readxl)
library(dplyr)
library(caret)
library(xgboost)
library(recipes)
library(themis)
library(tidyverse)
library(pROC)

# Set seed for reproducibility
set.seed(123)

# Set working directory
setwd("E:/CAG/Variant_Disease_Prediction_AI")

# Step 1: Load Data with error handling
gene_panel_file <- "HC_gene.xlsx"
clinvar_file <- "CV_HC.xlsx"

gene_panel <- tryCatch(
  read_excel(gene_panel_file),
  error = function(e) stop("Error reading gene panel file: ", e)
)

clinvar_data <- tryCatch(
  read_excel(clinvar_file),
  error = function(e) stop("Error reading ClinVar file: ", e)
)

# Step 2: Data Quality Checks and Filtering
# Print initial data summary
print("Initial data dimensions:")
print(dim(clinvar_data))
print("\nMissing values in ClinVar data:")
print(colSums(is.na(clinvar_data)))

# Filter ClinVar Data
genes_of_interest <- unique(gene_panel$Gene)
print(paste("Number of unique genes in panel:", length(genes_of_interest)))

filtered_clinvar <- clinvar_data %>%
  filter(Gene %in% genes_of_interest) %>%
  filter(!is.na(Gene), !is.na(Germline_classification))

# Step 3: Create Class variable and handle missing values
filtered_clinvar <- filtered_clinvar %>%
  mutate(Class = case_when(
    Germline_classification %in% c("Pathogenic", "Likely Pathogenic") ~ 1,
    Germline_classification %in% c("Benign", "Likely Benign") ~ 0,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(Class)) %>%
  mutate(Class = as.numeric(Class))  # Keep as numeric for XGBoost

print("\nClass distribution:")
print(table(filtered_clinvar$Class))

# Step 4: Create preprocessing recipe
prep_recipe <- recipe(Class ~ ., data = filtered_clinvar) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_rm(Germline_classification) %>%
  step_impute_mean(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_zv(all_predictors())

# Prepare the data
prepared_data <- prep(prep_recipe)
baked_data <- bake(prepared_data, new_data = NULL)

# Step 5: Split Data
train_index <- createDataPartition(baked_data$Class, p = 0.7, list = FALSE)
train_data <- baked_data[train_index, ]
test_data <- baked_data[-train_index, ]

# Store feature names for later use
train_features <- colnames(train_data %>% select(-Class))

# Step 6: Train XGBoost Model
# Create DMatrix objects
dtrain <- xgb.DMatrix(
  data = as.matrix(train_data %>% select(-Class)),
  label = train_data$Class
)

dtest <- xgb.DMatrix(
  data = as.matrix(test_data %>% select(-Class)),
  label = test_data$Class
)

# Set XGBoost parameters
xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = c("auc", "logloss"),
  max_depth = 6,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  nthread = parallel::detectCores() - 1
)

# Train with cross-validation
xgb_cv <- xgb.cv(
  params = xgb_params,
  data = dtrain,
  nrounds = 1000,
  nfold = 5,
  early_stopping_rounds = 50,
  verbose = TRUE
)

# Get optimal number of rounds
final_nrounds <- xgb_cv$best_iteration

# Train final model
xgb_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = final_nrounds,
  watchlist = list(train = dtrain, test = dtest),
  verbose = 1
)

# Step 7: Make Predictions
xgb_preds <- predict(xgb_model, dtest)

# Step 8: Evaluate Model
evaluate_model <- function(predictions, actual, threshold = 0.5) {
  pred_class <- as.factor(ifelse(predictions > threshold, 1, 0))
  actual <- as.factor(actual)
  
  # Ensure both have the same levels
  levels(pred_class) <- levels(actual) <- c("0", "1")
  
  # Create confusion matrix
  confusion_mat <- confusionMatrix(pred_class, actual)
  
  # Calculate metrics
  precision <- confusion_mat$byClass["Pos Pred Value"]
  recall <- confusion_mat$byClass["Sensitivity"]
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  # Calculate AUC
  auc <- roc(actual, predictions)$auc
  
  return(list(
    confusion_matrix = confusion_mat,
    accuracy = confusion_mat$overall["Accuracy"],
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    auc = auc
  ))
}

# Get evaluation metrics
model_evaluation <- evaluate_model(xgb_preds, test_data$Class)

# Print results
print("\nModel Evaluation Results:")
print("-------------------------")
print(paste("Accuracy:", round(model_evaluation$accuracy, 4)))
print(paste("Precision:", round(model_evaluation$precision, 4)))
print(paste("Recall:", round(model_evaluation$recall, 4)))
print(paste("F1 Score:", round(model_evaluation$f1_score, 4)))
print(paste("AUC:", round(model_evaluation$auc, 4)))
print("\nConfusion Matrix:")
print(model_evaluation$confusion_matrix$table)

# Step 9: Feature Importance
importance_matrix <- xgb.importance(
  feature_names = train_features,
  model = xgb_model
)

print("\nTop 10 Important Features:")
print(head(importance_matrix, ))

# Step 10: Save Results
# Save model
saveRDS(xgb_model, "xgboost_model.rds")
saveRDS(prepared_data, "preprocessing_recipe.rds")  # Save preprocessing recipe
write.csv(importance_matrix, "feature_importance.csv")

# Create results dataframe
results_df <- data.frame(
  actual = test_data$Class,
  predicted_prob = xgb_preds,
  predicted_class = ifelse(xgb_preds > 0.5, 1, 0)
)

write.csv(results_df, "prediction_results.csv", row.names = FALSE)

print("\nProcessing complete. Model and results have been saved.")
# Extract feature importance
importance_matrix <- xgb.importance(model = xgb_model, feature_names =  train_features)

# Plot feature importance
xgb.plot.importance(importance_matrix, top_n = 10, main = "Top 10 Feature Importances")
