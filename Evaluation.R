
# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(nlme)
library(glmertree)
library(lme4) 
library(partykit)
library(gpboost)
library(knitr)
library(REEMtree)
library(party)


# Load the data
all_data <- read.csv("D:/payanname/Matlabcodes/HbO_Behavior_fnirs.csv")

# Convert variables to factors as described in the paper or implied
all_data $n_back <- factor(all_data$n_back, levels = c("0-back", "2-back", "3-back"), labels = c("0", "2", "3")) 
all_data $Subject <- as.factor(all_data$Subject)
all_data $Indices <- as.factor(all_data$Indices)
all_data $Session <- as.factor(all_data$Session) 
all_data $Gender <- as.factor(all_data$Gender)  

# Filter for 2-back and 3-back conditions as per the paper
df_filtered <- all_data %>%
  filter(n_back %in% c("2", "3")) %>% # Use new labels if changed above
  droplevels() # Drop unused factor levels after filtering

# Define fixed and random effects based on the parsimonious model in the paper
# The paper states: "We have n-back and MeanRT variables in the parsimonious model."
fixed_effects_formula <- value ~ n_back + MeanRT
random_effects_formula <- ~ 1 | Subject/Indices

# Control parameters for lme, especially for potentially complex models
ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", msVerbose = FALSE)

# --- 1. LMM without nested structure ---
cat("Fitting LMM without nested structure...\n")
time_lmm_no_nest <- system.time({
  lmm_no_nest <- lme(fixed = fixed_effects_formula, # value ~ n_back + MeanRT
                     random = ~ 1 | Subject,       # Only Subject random effect
                     data = df_filtered,
                     method = "ML", 
                     control = lmeControl(opt = "optim")) 
})

actual_values <- df_filtered$value
predicted_lmm_no_nest <- fitted(lmm_no_nest)

mse_lmm_no_nest  <- mean((actual_values - predicted_lmm_no_nest)^2)
rmse_lmm_no_nest <- sqrt(mse_lmm_no_nest)
mae_lmm_no_nest  <- mean(abs(actual_values - predicted_lmm_no_nest))
time_val_lmm_no_nest <- time_lmm_no_nest[["elapsed"]] 

cat("Done.\n")
cat("Time taken for LMM without nested structure:", round(time_val_lmm_no_nest, 2), "s\n")



cat("Re-fitting LMM with nested structure (m4 with corSymm) to get accurate time...\n")
ctrl_m4 <- lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim", niterEM = 100, msVerbose = FALSE)

time_m4 <- system.time({
  m4_timed <- lme(fixed = fixed_effects_formula, # value ~ n_back + MeanRT
                  random = random_effects_formula, # ~1 | Subject/Indices
                  correlation = corSymm(form = ~ 1 | Subject/Indices),
                  data = df_filtered,
                  method = "ML",
                  control = ctrl_m4)
})

predicted_m4 <- fitted(m4_timed)

mse_m4  <- mean((actual_values - predicted_m4)^2)
rmse_m4 <- sqrt(mse_m4)
mae_m4  <- mean(abs(actual_values - predicted_m4))
time_val_m4 <- time_m4[["elapsed"]]

cat("Done.\n")
cat("Time taken for LMM with nested structure (corSymm):", round(time_val_m4, 2), "s\n")


results_table <- data.frame(
  Method = c("LMM without nested structure",
             "LMM with nested structure (corSymm)"),
  RMSE = c(rmse_lmm_no_nest, rmse_m4),
  MSE = c(mse_lmm_no_nest, mse_m4),
  MAE = c(mae_lmm_no_nest, mae_m4),
  System_time = c(time_val_lmm_no_nest, time_val_m4)
)

cat("\n--- Calculated Performance Metrics Table ---\n")
print(results_table, row.names = FALSE)


##########################################################################
#########2:GLMM tree

file_path <- "D:\\payanname\\Matlabcodes\\HbO_Behavior_fnirs.csv"
all_data <- read_csv(file_path, show_col_types = FALSE)

# Filter for 2-back and 3-back
data_for_glmm <- all_data %>%
  filter(n_back %in% c("2-back", "3-back")) %>%
  mutate(
    Subject = as.factor(Subject),
    Session = as.factor(Session),
    Gender = as.factor(Gender),
    n_back = factor(n_back, levels = c("2-back", "3-back")),
    Indices = as.factor(Indices)
  )

### STEP 2: FITTING GLMM TREES
# Base LMM tree model
lt <- lmertree(value ~ n_back + MeanRT + Session + Accuracy | Subject | Age + Gender,
               data = data_for_glmm)


# For comparison (as in article): define two variations of models
glmm_tree_subj_only <- lt  # current model without nesting

# Nested version (simulate for now using Session nested in Subject)
glmm_tree_nested <- lmertree(value ~ n_back + MeanRT + Session + Accuracy | (1 | Subject/Indices)  | Age + Gender,
                             data = data_for_glmm)



### STEP 3: MODEL EVALUATION
cat("--- Step 3: Evaluating Models ---\n")

calculate_metrics <- function(model, data, response_var_name = "value") {
  predictions <- predict(model, newdata = data, type = "response", allow.new.levels = TRUE)
  actual <- data[[response_var_name]]
  valid_indices <- !is.na(predictions) & !is.na(actual)
  mse <- mean((actual[valid_indices] - predictions[valid_indices])^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(actual[valid_indices] - predictions[valid_indices]))
  return(list(MSE = mse, RMSE = rmse, MAE = mae))
}

metrics_glmm_subj_only <- calculate_metrics(glmm_tree_subj_only, data_for_glmm)
metrics_glmm_nested <- calculate_metrics(glmm_tree_nested, data_for_glmm)

print("Metrics for model without nested structure:"); print(metrics_glmm_subj_only)
print("Metrics for model with nested structure:"); print(metrics_glmm_nested)


####################################
#######3:REEM tree

data_filtered_tbl <- all_data %>%
  filter(n_back %in% c("2-back", "3-back"))

data_filtered_tbl$Subject <- as.factor(data_filtered_tbl$Subject)
data_filtered_tbl$Session <- as.factor(data_filtered_tbl$Session)
data_filtered_tbl$Gender <- as.factor(data_filtered_tbl$Gender)
data_filtered_tbl$n_back <- as.factor(data_filtered_tbl$n_back)
data_filtered_tbl$Indices <- as.factor(data_filtered_tbl$Indices)

data_filtered_df <- as.data.frame(data_filtered_tbl)


fixed_effects_formula <- value ~ MeanRT + n_back + Age + Accuracy + Session + Gender
random_effects_formula <- ~ 1 | Subject


results_table <- data.frame(
  Structure = character(),
  Correlation = character(),
  RMSE = numeric(),
  MSE = numeric(),
  MAE = numeric(),
  stringsAsFactors = FALSE
)

compute_metrics <- function(model, newdata, actual) {
  preds <- predict(model, newdata = newdata)
  residuals <- actual - preds
  MSE <- mean(residuals^2, na.rm = TRUE)
  RMSE <- sqrt(MSE)
  MAE <- mean(abs(residuals), na.rm = TRUE)
  return(c(RMSE = RMSE, MSE = MSE, MAE = MAE))
}

# -----------------------
# Without nested - Default
model1 <- REEMtree(
  formula = value ~ MeanRT + n_back + Age + Accuracy + Session + Gender,
  data = data_filtered_df,
  random = ~1 | Subject
)
m1 <- compute_metrics(model1, data_filtered_df, data_filtered_df$value)
results_table <- rbind(results_table, data.frame(Structure = "Without nested", Correlation = "Default", t(m1)))

# -----------------------
# With nested - Default
data_filtered_df$Subject_Indices <- interaction(data_filtered_df$Subject, data_filtered_df$Indices)

model2 <- REEMtree(
  formula = value ~ MeanRT + n_back + Age + Accuracy + Session + Gender,
  data = data_filtered_df,
  random = ~1 | Subject_Indices
)
m2 <- compute_metrics(model2, data_filtered_df, data_filtered_df$value)
results_table <- rbind(results_table, data.frame(Structure = "With nested", Correlation = "Default", t(m2)))


print_scientific_table <- function(df, digits = 2) {
  df_formatted <- df
  df_formatted$RMSE <- formatC(df$RMSE, format = "e", digits = digits)
  df_formatted$MSE  <- formatC(df$MSE,  format = "e", digits = digits)
  df_formatted$MAE  <- formatC(df$MAE,  format = "e", digits = digits)
  print(df_formatted, row.names = FALSE)
}

print_scientific_table(results_table, digits = 2)

############################################
####################4:unbiased REEMtree

#-------------------------start of main function----------------------------------------------------
REEMctree_function <- function (formula, data, random, subset = NULL, initialRandomEffects = rep(0, 
                                                                                                 TotalObs), ErrorTolerance = 0.001, MaxIterations = 1000, 
                                verbose = FALSE, lme.control = lmeControl(returnObject = TRUE), ctree.control = party::ctree_control(mincriterion = 0.85),
                                method = "REML", correlation = NULL) 
{
  TotalObs <- dim(data)[1]
  if (identical(subset, NULL)) {
    subs <- rep(TRUE, dim(data)[1])
  }
  else {
    subs <- subset
  }
  Predictors <- paste(attr(terms(formula), "term.labels"), 
                      collapse = "+")
  TargetName <- formula[[2]]
  if (length(TargetName) > 1) 
    TargetName <- TargetName[3]
  if (verbose) 
    print(paste("Target variable: ", TargetName))
  Target <- data[, toString(TargetName)] 
  ContinueCondition <- TRUE
  iterations <- 0
  AdjustedTarget <- Target - initialRandomEffects
  oldlik <- -Inf
  newdata <- data
  newdata[, "SubsetVector"] <- subs
  while (ContinueCondition) {
    newdata[, "AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations + 1
    
    tree <- party::ctree(formula(paste(c("AdjustedTarget", Predictors), collapse = "~")), data = newdata, subset = subs, controls = ctree.control)  
    
    if (verbose) 
      print(tree)
    newdata[, "nodeInd"] <- 0 
    
    if(sum(subs) > 0 && !is.null(tree) && length(tree@tree$nodeID) > 0) { 
      newdata[subs, "nodeInd"] <- party::where(tree) 
    } else if (sum(subs) > 0) { 
      newdata[subs, "nodeInd"] <- 1 
      warning("ctree did not produce a valid tree structure in an iteration. All subset data assigned to node 1.")
    }
    
    if (any(is.na(newdata[subs, "nodeInd"]))) {
      stop("NA values in nodeInd after party::where call for subsetted data.")
    }
    
    current_node_indices_in_subset <- newdata[subs, "nodeInd"]
    
    if (min(current_node_indices_in_subset) == max(current_node_indices_in_subset)) { #it doesn't split on root for the current subset
      lmefit <- nlme::lme(formula(paste(c(toString(TargetName), 
                                          1), collapse = "~")), data = newdata, random = random, 
                          subset = SubsetVector, method = method, control = lme.control, 
                          correlation = correlation)
    }
    else {
      lmefit <- nlme::lme(formula(paste(c(toString(TargetName), 
                                          "as.factor(nodeInd)"), collapse = "~")), data = newdata, 
                          random = random, subset = SubsetVector, method = method,
                          control = lme.control, correlation = correlation)
    }
    if (verbose) {
      print(lmefit)
      print(paste("Estimated Error Variance = ", lmefit$sigma))
      print("Estimated Random Effects Variance = ")
      if(!is.null(lmefit$modelStruct$reStruct[[1]])) {
        print(as.matrix(lmefit$modelStruct$reStruct[[1]]) * 
                lmefit$sigma^2)
      } else {
        print("No random effects structure found or variance components are zero.")
      }
    }
    newlik <- logLik(lmefit)
    if (verbose) 
      print(paste("Log likelihood: ", newlik))
    
    if (is.na(newlik) || !is.finite(newlik)) {
      warning("Log-likelihood is NA or non-finite. Stopping iterations.")
      ContinueCondition <- FALSE
    } else {
      ContinueCondition <- (newlik - oldlik > ErrorTolerance & 
                              iterations < MaxIterations)
      oldlik <- newlik
    }
    
    AllEffects <- lmefit$residuals[, 1] - lmefit$residuals[, 
                                                           dim(lmefit$residuals)[2]]
    
    if(length(AllEffects) != sum(subs)) {
      stop("Mismatch in length of AllEffects and subsetted data. Check lme fit and residuals structure.")
    }
    AdjustedTarget[subs] <- Target[subs] - AllEffects
  }
  residuals <- rep(NA, length = length(Target))
  fitted_values_on_subset <- predict(lmefit) 
  if(length(fitted_values_on_subset) != sum(subs)) {
    stop("Mismatch in length of predictions and subsetted data.")
  }
  residuals[subs] <- Target[subs] - fitted_values_on_subset
  attr(residuals, "label") <- NULL 
  
  result <- list(Tree = tree, EffectModel = lmefit, RandomEffects = ranef(lmefit), 
                 BetweenMatrix = if(!is.null(lmefit$modelStruct$reStruct[[1]])) {
                   as.matrix(lmefit$modelStruct$reStruct[[1]]) * lmefit$sigma^2
                 } else { NA }, 
                 ErrorVariance = lmefit$sigma^2, data = data, 
                 logLik = newlik, IterationsUsed = iterations, Formula = formula, 
                 Random = random, Subset = subs, ErrorTolerance = ErrorTolerance, 
                 correlation = correlation, residuals = residuals, method = method, 
                 lme.control = lme.control, ctree.control = ctree.control)
  class(result) <- "REEMctree"
  return(result)
}
#-------------------------end of main function----------------------------------------------------
#--------------------------------------------------------------------------------------------------

##################################################


data_filtered_tbl <- all_data %>%
  filter(n_back %in% c("2-back", "3-back"))

data_filtered_tbl$Subject <- as.factor(data_filtered_tbl$Subject)
data_filtered_tbl$Session <- as.factor(data_filtered_tbl$Session)
data_filtered_tbl$Gender <- as.factor(data_filtered_tbl$Gender)
data_filtered_tbl$n_back <- as.factor(data_filtered_tbl$n_back)
data_filtered_tbl$Indices <- as.factor(data_filtered_tbl$Indices)

data_filtered_df <- as.data.frame(data_filtered_tbl)

fixed_effects_formula <- value ~ MeanRT + n_back + Age + Accuracy + Session + Gender
data_filtered_df$Subject_Indices <- interaction(data_filtered_df$Subject, data_filtered_df$Indices, drop = TRUE)


print("--- Starting REEMctree model: Default (Without nested) ---")
start_time_ctree_simple <- Sys.time()
reemctree_model_simple <- tryCatch({
  REEMctree_function(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = ~1 | Subject,
    ctree.control = party::ctree_control(mincriterion = 0.85)
  )
}, error = function(e) { NULL })
system_time_ctree_simple <- Sys.time() - start_time_ctree_simple


metrics_ctree_simple <- list(MSE = NA, RMSE = NA, MAE = NA)
if (!is.null(reemctree_model_simple)) {
  residuals_from_model <- reemctree_model_simple$residuals[reemctree_model_simple$Subset]
  metrics_ctree_simple$MSE <- mean(residuals_from_model^2, na.rm = TRUE)
  metrics_ctree_simple$RMSE <- sqrt(metrics_ctree_simple$MSE)
  metrics_ctree_simple$MAE <- mean(abs(residuals_from_model), na.rm = TRUE)
}


print("--- Starting REEMctree model: Default (With nested) ---")
start_time_ctree_nested <- Sys.time()
reemctree_model_nested <- tryCatch({
  REEMctree_function(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = ~1 | Subject_Indices,
    ctree.control = party::ctree_control(mincriterion = 0.85)
  )
}, error = function(e) { NULL })
system_time_ctree_nested <- Sys.time() - start_time_ctree_nested


metrics_ctree_nested <- list(MSE = NA, RMSE = NA, MAE = NA)
if (!is.null(reemctree_model_nested)) {
  residuals_from_model_nested <- reemctree_model_nested$residuals[reemctree_model_nested$Subset]
  metrics_ctree_nested$MSE <- mean(residuals_from_model_nested^2, na.rm = TRUE)
  metrics_ctree_nested$RMSE <- sqrt(metrics_ctree_nested$MSE)
  metrics_ctree_nested$MAE <- mean(abs(residuals_from_model_nested), na.rm = TRUE)
}


final_default_table <- data.frame(
  Method_Structure = character(),
  Correlation_Type = character(),
  RMSE = numeric(),
  MSE = numeric(),
  MAE = numeric(),
  System_Time_sec = numeric(),
  stringsAsFactors = FALSE
)


if(!is.na(metrics_ctree_simple$RMSE)) {
  final_default_table <- rbind(final_default_table, data.frame(
    Method_Structure = "Unbiased RE-EM tree - Without nested structure",
    Correlation_Type = "Default",
    RMSE = metrics_ctree_simple$RMSE,
    MSE = metrics_ctree_simple$MSE,
    MAE = metrics_ctree_simple$MAE,
    System_Time_sec = as.numeric(system_time_ctree_simple, units = "secs")
  ))
}


if(!is.na(metrics_ctree_nested$RMSE)) {
  final_default_table <- rbind(final_default_table, data.frame(
    Method_Structure = "Unbiased RE-EM tree - With nested structure",
    Correlation_Type = "Default",
    RMSE = metrics_ctree_nested$RMSE,
    MSE = metrics_ctree_nested$MSE,
    MAE = metrics_ctree_nested$MAE,
    System_Time_sec = as.numeric(system_time_ctree_nested, units = "secs")
  ))
}

if (nrow(final_default_table) > 0) {
  formatted_table <- final_default_table
  formatted_table$RMSE <- sprintf("%.1e", formatted_table$RMSE)
  formatted_table$MSE <- sprintf("%.1e", formatted_table$MSE)
  formatted_table$MAE <- sprintf("%.1e", formatted_table$MAE)
  formatted_table$System_Time_sec <- sprintf("%.2f", formatted_table$System_Time_sec)
  
  cat("\n\n--- Performance metrics for Unbiased RE-EM tree (Default Models Only) ---\n")
  print(formatted_table, row.names = FALSE)
} else {
  cat("\n\n--- No valid model results to populate the performance metrics table. ---\n")
}
####################################################
#################5:gpboost

all_data <- read.csv("D:/payanname/Matlabcodes/HbO_Behavior_fnirs.csv")

data <- all_data %>%
  filter(n_back %in% c("2-back", "3-back")) %>%
  mutate(
    Subject = as.factor(Subject), Session = as.factor(Session),
    Gender = as.factor(Gender), n_back = factor(n_back, levels = c("2-back", "3-back")),
    Indices = as.factor(Indices)
  ) %>%
  mutate(
    Age = as.numeric(scale(Age)),
    Accuracy = as.numeric(scale(Accuracy)),
    MeanRT = as.numeric(scale(MeanRT))
  )

# STEP 2: PREPARE DATA FOR GPBOOST 
X <- model.matrix(value ~ Age + Accuracy + MeanRT + n_back + Session + Gender - 1, data = data)
y <- data$value
SCALING_FACTOR <- 1e5
y_scaled <- y * SCALING_FACTOR
dtrain <- gpb.Dataset(data = X, label = y_scaled)
results_table4 <- data.frame()

# ===================================================================
# MODEL 1: GPBoost without nested structure (using WLS/newton)
# ===================================================================
cat("Running GPBoost WITHOUT nested structure using WLS (newton) optimizer...\n")

group_data_wn <- data$Subject
gp_model_wn <- GPModel(group_data = group_data_wn, likelihood = "gaussian")
# Explicitly set the optimizer to 'newton' which corresponds to WLS
gp_model_wn$set_optim_params(params = list(optimizer_cov = "newton"))

start_time <- proc.time()
bst_wn <- gpboost(
  data = dtrain, gp_model = gp_model_wn,
  nrounds = 4, learning_rate = 0.1, verbose = 0
)
end_time <- proc.time()
exec_time <- (end_time - start_time)["user.self"]

preds_scaled <- predict(bst_wn, data = X, group_data_pred = group_data_wn)$response_mean
preds_unscaled <- preds_scaled / SCALING_FACTOR

mse <- mean((y - preds_unscaled)^2)
rmse <- sqrt(mse)
mae <- mean(abs(y - preds_unscaled))

results_table4 <- rbind(results_table4, data.frame(
  Method = "GPBoost without nested structure",
  RMSE = rmse, MSE = mse, MAE = mae, System_Time = exec_time
))

# ===================================================================
# MODEL 2: GPBoost with nested structure (using WLS/newton)
# ===================================================================
cat("Running GPBoost WITH nested structure using WLS (newton) optimizer...\n")

group_data_n <- as.matrix(data[, c("Subject", "Indices")])
gp_model_n <- GPModel(group_data = group_data_n, likelihood = "gaussian")
# Explicitly set the optimizer to 'newton'
gp_model_n$set_optim_params(params = list(optimizer_cov = "newton"))

start_time <- proc.time()
bst_n <- gpboost(
  data = dtrain, gp_model = gp_model_n,
  nrounds = 6, learning_rate = 0.1, verbose = 0
)
end_time <- proc.time()
exec_time <- (end_time - start_time)["user.self"]

preds_scaled_n <- predict(bst_n, data = X, group_data_pred = group_data_n)$response_mean
preds_unscaled_n <- preds_scaled_n / SCALING_FACTOR

mse_n <- mean((y - preds_unscaled_n)^2)
rmse_n <- sqrt(mse_n)
mae_n <- mean(abs(y - preds_unscaled_n))

results_table4 <- rbind(results_table4, data.frame(
  Method = "GPBoost with nested structure",
  RMSE = rmse_n, MSE = mse_n, MAE = mae_n, System_Time = exec_time
))

# STEP 4: DISPLAY THE FINAL TABLE (WITHOUT WARNINGS)
formatted_results <- results_table4 %>%
  mutate(
    RMSE = format(RMSE, scientific = TRUE, digits = 2),
    MSE = format(MSE, scientific = TRUE, digits = 2),
    MAE = format(MAE, scientific = TRUE, digits = 2),
    System_Time = round(System_Time, 2)
  )

# Use suppressWarnings to hide the harmless kable/xfun warnings
suppressWarnings({
  print(kable(formatted_results, caption = "Replication of GPBoost Rows from Table 4 (with WLS optimizer)"))
})


#############################################

# ==========
# Table 6 
# ==========


REEMctree_function <- function (formula, data, random, subset = NULL, initialRandomEffects = rep(0, nrow(data)), ErrorTolerance = 0.001, MaxIterations = 1000, 
                                verbose = FALSE, lme.control = lmeControl(returnObject = TRUE), ctree.control = party::ctree_control(mincriterion = 0.85),
                                method = "REML", correlation = NULL) 
{
  if (is.null(subset)) {
    subs <- rep(TRUE, nrow(data))
  } else {
    subs <- subset
  }
  
  Predictors <- paste(attr(terms(formula), "term.labels"), collapse = "+")
  TargetName <- deparse(formula[[2]])
  Target <- data[[TargetName]]
  
  ContinueCondition <- TRUE
  iterations <- 0
  AdjustedTarget <- Target - initialRandomEffects
  oldlik <- -Inf
  
  newdata <- data
  # *** FIX IS HERE: The column for subsetting MUST be created ***
  newdata$SubsetVector <- subs
  
  while (ContinueCondition) {
    newdata$AdjustedTarget <- AdjustedTarget
    iterations <- iterations + 1
    
    tree <- party::ctree(formula(paste("AdjustedTarget ~", Predictors)), data = newdata, subset = subs, controls = ctree.control)
    
    newdata$nodeInd <- as.factor(1)
    if(sum(subs) > 0 && length(tree@tree) > 0) {
      newdata$nodeInd[subs] <- as.factor(party::where(tree))
    }
    
    if (length(unique(newdata$nodeInd[subs])) == 1) {
      lmefit <- lme(formula(paste(TargetName, "~ 1")), data = newdata, random = random, 
                    subset = SubsetVector, method = method, control = lme.control, 
                    correlation = correlation)
    } else {
      lmefit <- lme(formula(paste(TargetName, "~ nodeInd")), data = newdata, 
                    random = random, subset = SubsetVector, method = method, 
                    control = lme.control, correlation = correlation)
    }
    
    newlik <- logLik(lmefit)
    
    if (is.na(newlik) || !is.finite(newlik) || (newlik - oldlik < ErrorTolerance) || (iterations >= MaxIterations)) {
      ContinueCondition <- FALSE
    } else {
      oldlik <- newlik
      AllEffects <- fitted(lmefit, level = 0) - fitted(lmefit, level = length(lmefit$groups))
      AdjustedTarget[subs] <- Target[subs] - AllEffects
    }
  }
  
  residuals <- rep(NA, length(Target))
  residuals[subs] <- Target[subs] - predict(lmefit)
  
  return(list(residuals = residuals))
}

all_data <- read.csv("D:/payanname/Matlabcodes/HbO_Behavior_fnirs.csv")
df_original <- all_data %>%
  filter(n_back %in% c("2-back", "3-back")) %>%
  mutate(across(c(Subject, Session, Gender, n_back, Indices), as.factor)) %>%
  droplevels()
df_original$Subject_Indices <- interaction(df_original$Subject, df_original$Indices, drop = TRUE)

set.seed(12)
sd_value_original <- sd(df_original$value, na.rm = TRUE)
n_to_contaminate <- round(0.10 * nrow(df_original))
indices_to_contaminate <- sample(1:nrow(df_original), n_to_contaminate)

df_contaminated_10sd <- df_original
df_contaminated_10sd$value[indices_to_contaminate] <- df_contaminated_10sd$value[indices_to_contaminate] + (10 * sd_value_original)
df_contaminated_10sd$Subject_Indices <- interaction(df_contaminated_10sd$Subject, df_contaminated_10sd$Indices, drop = TRUE)

df_contaminated_15sd <- df_original
df_contaminated_15sd$value[indices_to_contaminate] <- df_contaminated_15sd$value[indices_to_contaminate] + (15 * sd_value_original)
df_contaminated_15sd$Subject_Indices <- interaction(df_contaminated_15sd$Subject, df_contaminated_15sd$Indices, drop = TRUE)


final_table6_results <- data.frame()
calc_rmse <- function(res) { if(is.null(res) || all(is.na(res))) NA else sqrt(mean(res^2, na.rm=TRUE)) }
calc_mae <- function(res) { if(is.null(res) || all(is.na(res))) NA else mean(abs(res), na.rm=TRUE) }

model_list <- c("LMM", "GLMMtree", "REEMtree", "Unbiased REEMtree", "GPBoost")

for (model_name in model_list) {
  cat(paste("\nProcessing Model:", model_name, "...\n"))
  tryCatch({
    get_residuals <- function(model_name, df) {
      suppressWarnings({
        if(model_name == "LMM") {
          nn <- lme(value ~ n_back + MeanRT, random = ~1|Subject, data = df)
          n <- lme(value ~ n_back + MeanRT, random = ~1|Subject/Indices, data = df)
          return(list(nn = residuals(nn), n = residuals(n)))
        } else if (model_name == "GLMMtree") {
          nn <- lmertree(value ~ n_back + MeanRT | Subject | Age + Gender, data = df)
          n  <- lmertree(value ~ n_back + MeanRT | (1|Subject/Indices) | Age + Gender, data = df)
          return(list(nn = df$value - predict(nn), n = df$value - predict(n)))
        } else if (model_name == "REEMtree") {
          nn <- REEMtree(value ~ MeanRT + n_back, data=df, random=~1|Subject)
          n  <- REEMtree(value ~ MeanRT + n_back, data=df, random=~1|Subject_Indices)
          return(list(nn = nn$residuals, n = n$residuals))
        } else if (model_name == "Unbiased REEMtree") {
          nn <- REEMctree_function(value ~ MeanRT + n_back, data=df, random=~1|Subject)
          n  <- REEMctree_function(value ~ MeanRT + n_back, data=df, random=~1|Subject_Indices)
          return(list(nn = nn$residuals, n = n$residuals))
        } else if (model_name == "GPBoost") {
          X <- model.matrix(value~Age+Accuracy+MeanRT+n_back+Session+Gender-1, df)
          y_scaled <- df$value * 1e5
          d <- gpb.Dataset(X, label=y_scaled)
          gpm_nn <- GPModel(group_data=df$Subject); bst_nn <- gpboost(data=d, gp_model=gpm_nn, nrounds=4, verbose=0)
          gpm_n <- GPModel(group_data=as.matrix(df[,c("Subject","Indices")])); bst_n <- gpboost(data=d, gp_model=gpm_n, nrounds=6, verbose=0)
          pred_nn <- predict(bst_nn, data=X, group_data_pred=df$Subject)$response_mean/1e5
          pred_n <- predict(bst_n, data=X, group_data_pred=as.matrix(df[,c("Subject","Indices")]))$response_mean/1e5
          return(list(nn = df$value - pred_nn, n = df$value - pred_n))
        }
      })
    }
    
    res_orig <- get_residuals(model_name, df_original)
    res_c10  <- get_residuals(model_name, df_contaminated_10sd)
    res_c15  <- get_residuals(model_name, df_contaminated_15sd)
    
    rmse_orig_nn <- calc_rmse(res_orig$nn); mae_orig_nn <- calc_mae(res_orig$nn)
    rmse_c10_nn <- calc_rmse(res_c10$nn); mae_c10_nn <- calc_mae(res_c10$nn)
    rmse_c15_nn <- calc_rmse(res_c15$nn); mae_c15_nn <- calc_mae(res_c15$nn)
    final_table6_results <- rbind(final_table6_results, data.frame(
      Methods = paste(model_name, "without nested structure"),
      RMSE_Cont1 = ((rmse_c10_nn - rmse_orig_nn) / rmse_orig_nn) * 100,
      MAE_Cont1 = ((mae_c10_nn - mae_orig_nn) / mae_orig_nn) * 100,
      RMSE_Cont2 = ((rmse_c15_nn - rmse_orig_nn) / rmse_orig_nn) * 100,
      MAE_Cont2 = ((mae_c15_nn - mae_orig_nn) / mae_orig_nn) * 100))
    
    rmse_orig_n <- calc_rmse(res_orig$n); mae_orig_n <- calc_mae(res_orig$n)
    rmse_c10_n <- calc_rmse(res_c10$n); mae_c10_n <- calc_mae(res_c10$n)
    rmse_c15_n <- calc_rmse(res_c15$n); mae_c15_n <- calc_mae(res_c15$n)
    final_table6_results <- rbind(final_table6_results, data.frame(
      Methods = paste(model_name, "with nested structure"),
      RMSE_Cont1 = ((rmse_c10_n - rmse_orig_n) / rmse_orig_n) * 100,
      MAE_Cont1 = ((mae_c10_n - mae_orig_n) / mae_orig_n) * 100,
      RMSE_Cont2 = ((rmse_c15_n - rmse_orig_n) / rmse_orig_n) * 100,
      MAE_Cont2 = ((mae_c15_n - mae_orig_n) / mae_orig_n) * 100))
    
  }, error = function(e) {
    cat(paste("  -> ERROR processing", model_name, ":", e$message, "\n"))
  })
}

table6_formatted <- final_table6_results
table6_formatted[, -1] <- round(table6_formatted[, -1], 1)
colnames(table6_formatted) <- c("Methods", "RMSE (%) Contamination 1", "MAE (%) Contamination 1", "RMSE (%) Contamination 2", "MAE (%) Contamination 2")

suppressWarnings({
  print(knitr::kable(table6_formatted, caption = "Table 6: Relative Change (%) in Performance Metrics", align = 'lrrrr'))
})



######################################################cross validation
cat("\n\n--- Starting Cross-Validation Analysis (Replicating Table 7) ---\n")
set.seed(40)

# ===================================================================
# STEP 1: SCALING and DATA SPLIT 
# ===================================================================
SCALING_FACTOR <- 10000 
df_filtered_scaled <- df_filtered %>%
  mutate(value_scaled = value * SCALING_FACTOR)

unique_subjects <- unique(df_filtered_scaled$Subject)
train_subjects <- sample(unique_subjects, size = round(0.9 * length(unique_subjects))) # use round for safety

train_data <- df_filtered_scaled[df_filtered_scaled$Subject %in% train_subjects, ]
test_data  <- df_filtered_scaled[!(df_filtered_scaled$Subject %in% train_subjects), ]


cv_results_table <- data.frame()


calculate_test_metrics <- function(actual_values, predicted_values) {
  mse <- mean((actual_values - predicted_values)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  mae <- mean(abs(actual_values - predicted_values), na.rm = TRUE)
  return(list(RMSE = rmse, MSE = mse, MAE = mae))
}

# ===================================================================
# STEP 2: MODEL FITTING and EVALUATION
# ===================================================================

full_fixed_formula <- value_scaled ~ n_back + MeanRT + Session + Accuracy

# -----------------
# 1. LMM 
# -----------------
cat("\nProcessing LMM for Cross-Validation...\n")
library(lme4) 

# LMM without nested
lmm_cv_no_nest <- lmer(update(full_fixed_formula, . ~ . + (1 | Subject)), data = train_data)
preds_lmm_no_nest <- predict(lmm_cv_no_nest, newdata = test_data, re.form = NA) / SCALING_FACTOR
metrics <- calculate_test_metrics(test_data$value, preds_lmm_no_nest)
cv_results_table <- rbind(cv_results_table, data.frame(Method = "LMM without nested structure", t(unlist(metrics))))

# LMM with nested
lmm_cv_nest <- lmer(update(full_fixed_formula, . ~ . + (1 | Subject/Indices)), data = train_data)
preds_lmm_nest <- predict(lmm_cv_nest, newdata = test_data, re.form = NA) / SCALING_FACTOR
metrics <- calculate_test_metrics(test_data$value, preds_lmm_nest)
cv_results_table <- rbind(cv_results_table, data.frame(Method = "LMM with nested structure", t(unlist(metrics))))


# -----------------
# 2. GLMM Tree 
# -----------------
cat("Processing GLMMtree for Cross-Validation...\n")

# GLMMtree without nested
glmm_tree_cv_no_nest <- lmertree(
  value_scaled ~ n_back + MeanRT + Session + Accuracy | Subject | Age + Gender, 
  data = train_data, 
  maxdepth = 3 
)
preds_glmm_no_nest <- predict(glmm_tree_cv_no_nest, newdata = test_data, allow.new.levels = TRUE) / SCALING_FACTOR
metrics <- calculate_test_metrics(test_data$value, preds_glmm_no_nest)
cv_results_table <- rbind(cv_results_table, data.frame(Method = "GLMM tree without nested structure", t(unlist(metrics))))

# GLMMtree with nested
glmm_tree_cv_nest <- lmertree(
  value_scaled ~ n_back + MeanRT + Session + Accuracy | (1 | Subject/Indices) | Age + Gender, 
  data = train_data,
  maxdepth = 3 
)
preds_glmm_nest <- predict(glmm_tree_cv_nest, newdata = test_data, allow.new.levels = TRUE) / SCALING_FACTOR
metrics <- calculate_test_metrics(test_data$value, preds_glmm_nest)
cv_results_table <- rbind(cv_results_table, data.frame(Method = "GLMM tree with nested structure", t(unlist(metrics))))
# -----------------
# 3. RE-EM Tree
# -----------------
cat("Processing RE-EM tree for Cross-Validation...\n")

# REEMtree without nested
reem_cv_no_nest <- REEMtree(value ~ MeanRT + n_back + Age + Accuracy + Session + Gender, data = train_data, random = ~1 | Subject)

preds_reem_no_nest <- predict(reem_cv_no_nest$Tree, newdata = test_data) 
metrics <- calculate_test_metrics(test_data$value, preds_reem_no_nest)
cv_results_table <- rbind(cv_results_table, data.frame(Method = "RE-EM tree without nested structure", t(unlist(metrics))))

# REEMtree with nested
train_data$Subject_Indices <- interaction(train_data$Subject, train_data$Indices, drop = TRUE)
reem_cv_nest <- REEMtree(value ~ MeanRT + n_back + Age + Accuracy + Session + Gender, data = train_data, random = ~1 | Subject_Indices)
preds_reem_nest <- predict(reem_cv_nest$Tree, newdata = test_data)
metrics <- calculate_test_metrics(test_data$value, preds_reem_nest)
cv_results_table <- rbind(cv_results_table, data.frame(Method = "RE-EM tree with nested structure", t(unlist(metrics))))
# -----------------
# 4. GPBoost 
# -----------------
library(fastDummies)
cat("Processing GPBoost for Cross-Validation using fastDummies...\n")

cols_to_dummy <- c("n_back", "Session", "Gender")
numeric_cols <- c("Age", "Accuracy", "MeanRT")
all_cols_gp <- c(cols_to_dummy, numeric_cols)


combined_data_gp <- rbind(train_data[, all_cols_gp], test_data[, all_cols_gp])

combined_data_dummied_gp <- dummy_cols(combined_data_gp, select_columns = cols_to_dummy, remove_first_dummy = FALSE, remove_selected_columns = TRUE)

n_train <- nrow(train_data)
X_train <- as.matrix(combined_data_dummied_gp[1:n_train, ])
X_test  <- as.matrix(combined_data_dummied_gp[(n_train + 1):nrow(combined_data_dummied_gp), ])

y_train_scaled <- train_data$value * 1e5

# GPBoost without nested
gp_model_cv_wn <- GPModel(group_data=train_data$Subject, likelihood="gaussian")
bst_cv_wn <- gpboost(data=gpb.Dataset(X_train, label=y_train_scaled), gp_model=gp_model_cv_wn, nrounds=4, verbose=0)
dummy_group_pred_wn <- rep("new_subject", nrow(test_data))
preds_gp_wn <- predict(bst_cv_wn, data=X_test, group_data_pred=dummy_group_pred_wn)$response_mean / 1e5
cv_results_table <- rbind(cv_results_table, data.frame(Method="GPBoost without nested structure", t(unlist(calculate_test_metrics(test_data$value, preds_gp_wn)))))

# GPBoost with nested
group_data_train_n <- as.matrix(train_data[, c("Subject", "Indices")])
gp_model_cv_n <- GPModel(group_data=group_data_train_n, likelihood="gaussian")
bst_cv_n <- gpboost(data=gpb.Dataset(X_train, label=y_train_scaled), gp_model=gp_model_cv_n, nrounds=6, verbose=0)
dummy_group_pred_n <- cbind(Subject=rep("new_subject",nrow(test_data)), Indices=paste0("new_index_",1:nrow(test_data)))
preds_gp_n <- predict(bst_cv_n, data=X_test, group_data_pred=dummy_group_pred_n)$response_mean / 1e5
cv_results_table <- rbind(cv_results_table, data.frame(Method="GPBoost with nested structure", t(unlist(calculate_test_metrics(test_data$value, preds_gp_n)))))


cat("\n\n--- Final Cross-Validation Results (Table 7 Replication) ---\n")
formatted_cv_table <- cv_results_table
formatted_cv_table$RMSE <- formatC(cv_results_table$RMSE, format = "e", digits = 1)
formatted_cv_table$MSE <- formatC(cv_results_table$MSE, format = "e", digits = 1)
formatted_cv_table$MAE <- formatC(cv_results_table$MAE, format = "e", digits = 1)
print(formatted_cv_table, row.names = FALSE)






