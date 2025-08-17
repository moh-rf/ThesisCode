# Load following package first
library(party) 

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


library(readr)
library(dplyr)
library(nlme) 
library(ggplot2)
library(gridExtra)

all_data <- read_csv("D:/payanname/Matlabcodes/HbO_Behavior_fnirs.csv")

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




print("--- Starting REEMctree model: simple structure ---")
start_time_ctree_simple <- Sys.time()
reemctree_model_simple <- tryCatch({
  REEMctree_function(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = ~1 | Subject,
    ctree.control = party::ctree_control(mincriterion = 0.85)
  )
}, error = function(e) {
  print("Error running REEMctree (simple):")
  print(e)
  NULL
})
system_time_ctree_simple <- Sys.time() - start_time_ctree_simple

# Prediction and evaluation
metrics_ctree_simple <- list(MSE = NA, RMSE = NA, MAE = NA)
if (!is.null(reemctree_model_simple) && !is.null(reemctree_model_simple$EffectModel)) {

  actuals_ctree_simple <- data_filtered_df$value[reemctree_model_simple$Subset]
  predictions_ctree_simple <- actuals_ctree_simple - reemctree_model_simple$residuals[reemctree_model_simple$Subset]
  

  residuals_from_model <- reemctree_model_simple$residuals[reemctree_model_simple$Subset] 
  

  if(any(is.na(residuals_from_model))) {
    warning("NA values found in residuals_from_model (simple). Metrics might be affected.")
  }
  
  metrics_ctree_simple$MSE <- mean(residuals_from_model^2, na.rm = TRUE)
  metrics_ctree_simple$RMSE <- sqrt(metrics_ctree_simple$MSE)
  metrics_ctree_simple$MAE <- mean(abs(residuals_from_model), na.rm = TRUE)
} else if (!is.null(reemctree_model_simple) && is.null(reemctree_model_simple$EffectModel)) {
  print("REEMctree (simple) ran but EffectModel is NULL. Cannot calculate metrics.")
}


####################################################################################

print("--- Starting REEMctree model: nested structure ---")
start_time_ctree_nested <- Sys.time()
reemctree_model_nested <- tryCatch({
  REEMctree_function(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = ~1 | Subject_Indices,
    ctree.control = party::ctree_control(mincriterion = 0.85) 
  )
}, error = function(e) {
  print("Error running REEMctree (nested):")
  print(e)
  # traceback()
  NULL
})
system_time_ctree_nested <- Sys.time() - start_time_ctree_nested

# Prediction and evaluation
metrics_ctree_nested <- list(MSE = NA, RMSE = NA, MAE = NA)
if (!is.null(reemctree_model_nested) && !is.null(reemctree_model_nested$EffectModel)) {
  residuals_from_model_nested <- reemctree_model_nested$residuals[reemctree_model_nested$Subset]
  
  if(any(is.na(residuals_from_model_nested))) {
    warning("NA values found in residuals_from_model (nested). Metrics might be affected.")
  }
  
  metrics_ctree_nested$MSE <- mean(residuals_from_model_nested^2, na.rm = TRUE)
  metrics_ctree_nested$RMSE <- sqrt(metrics_ctree_nested$MSE)
  metrics_ctree_nested$MAE <- mean(abs(residuals_from_model_nested), na.rm = TRUE)
  

  actuals_ctree_nested_plot <- data_filtered_df$value[reemctree_model_nested$Subset]
  predictions_ctree_nested_plot <- actuals_ctree_nested_plot - residuals_from_model_nested
  
} else if (!is.null(reemctree_model_nested) && is.null(reemctree_model_nested$EffectModel)) {
  print("REEMctree (nested) ran but EffectModel is NULL. Cannot calculate metrics.")
}


cat("\n--- REEMctree without nested model results ---\n")
print(metrics_ctree_simple)
print(system_time_ctree_simple)

cat("\n--- REEMctree with nested model results ---\n") 
print(metrics_ctree_nested)
print(system_time_ctree_nested)

plot_ctree_simple <- ggplot() + theme_void()
if (!is.null(reemctree_model_simple) && !is.null(reemctree_model_simple$EffectModel) && !all(is.na(metrics_ctree_simple))) {
  actuals_ctree_simple_plot <- data_filtered_df$value[reemctree_model_simple$Subset]
  predictions_ctree_simple_plot <- actuals_ctree_simple_plot - reemctree_model_simple$residuals[reemctree_model_simple$Subset]
  df_simple <- data.frame(Fitted = predictions_ctree_simple_plot, Actual = actuals_ctree_simple_plot)
  
  plot_ctree_simple <- ggplot(df_simple, aes(x = Fitted, y = Actual)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(x = "Fitted Values", y = "Actual Values",
         caption = "(b) Unbiased RE-EM tree without nested data structure") +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 10, b = 5)),
      plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
    )
}

plot_ctree_nested <- ggplot() + theme_void()
if (!is.null(reemctree_model_nested) && !is.null(reemctree_model_nested$EffectModel) && !all(is.na(metrics_ctree_nested))) {
  df_nested <- data.frame(Fitted = predictions_ctree_nested_plot, Actual = actuals_ctree_nested_plot)
  
  plot_ctree_nested <- ggplot(df_nested, aes(x = Fitted, y = Actual)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(x = "Fitted Values", y = "Actual Values",
         caption = "(a) Unbiased RE-EM tree with nested data structure") +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 10, b = 5)),
      plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
    )
}

grid.arrange(plot_ctree_nested, plot_ctree_simple, ncol = 2)


##################################################

## different correlation structures for Unbiased RE-EM tree without nested structure
cor_structures <- list(
  "NoCorrelation" = NULL,
  "AR1" = nlme::corAR1(form = ~1 | Subject), 
  "CAR1" = nlme::corCAR1(form = ~1 | Subject), 
  "CompSymm" = nlme::corCompSymm(form = ~1 | Subject) 
)


get_loglik <- function(cor_struct_name, cor_structure_obj) { 
  print(paste("Running model for correlation:", cor_struct_name))
  model_output <- tryCatch({
    REEMctree_function(
      formula = fixed_effects_formula,
      data = data_filtered_df,
      random = ~1 | Subject,
      method = "ML", 
      correlation = cor_structure_obj, 
      ctree.control = party::ctree_control(mincriterion = 0.85)  
    )
  }, error = function(e) {
    print(paste("Error in get_loglik for", cor_struct_name, ":", e$message))
    return(NULL) 
  })
  
  if (!is.null(model_output) && !is.null(model_output$EffectModel)) {
    return(as.numeric(logLik(model_output$EffectModel)))
  } else {
    return(NA) 
  }
}


loglik_results <- sapply(names(cor_structures), function(name) get_loglik(name, cor_structures[[name]]))


loglik_df <- data.frame(
  Correlation_Structure = c("No within-group correlation", "AR(1)", "Continuous AR(1)", "Compound-symmetric"),
  Argument = c("Default", "corAR1()", "corCAR1()", "corCompSymm()"),
  LogLikelihood = round(as.numeric(loglik_results), 2) 
)

print(loglik_df)

#####################################################

## different correlation structures for Unbiased RE-EM tree with nested structure
cor_structures_nested <- list(
  "NoCorrelation" = NULL,
  "AR1" = nlme::corAR1(form = ~1 | Subject_Indices), 
  "CAR1" = nlme::corCAR1(form = ~1 | Subject_Indices), 
  "CompSymm" = nlme::corCompSymm(form = ~1 | Subject_Indices)
)


get_loglik_nested <- function(cor_struct_name, cor_structure_obj) {
  print(paste("Running nested model for correlation:", cor_struct_name))
  model_output <- tryCatch({
    REEMctree_function(
      formula = fixed_effects_formula,
      data = data_filtered_df,
      random = ~1 | Subject_Indices,
      method = "ML",
      correlation = cor_structure_obj,
      ctree.control = party::ctree_control(mincriterion = 0.85)  
    )
  }, error = function(e) {
    print(paste("Error in get_loglik_nested for", cor_struct_name, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(model_output) && !is.null(model_output$EffectModel)) {
    return(as.numeric(logLik(model_output$EffectModel)))
  } else {
    return(NA)
  }
}

loglik_results_nested <- sapply(names(cor_structures_nested), function(name) get_loglik_nested(name, cor_structures_nested[[name]]))

loglik_df_nested <- data.frame(
  Correlation_Structure = c("No within-group correlation", "AR(1)", "Continuous AR(1)", "Compound-symmetric"),
  Argument = c("Default", "corAR1()", "corCAR1()", "corCompSymm()"),
  LogLikelihood = round(as.numeric(loglik_results_nested), 2) 
)

print(loglik_df_nested)


##########################################


# Continuous AR(1)
print("--- Starting REEMctree model: Without Nested - Continuous AR(1) ---")
start_time_ctree_wn_car1 <- Sys.time()
reemctree_model_wn_car1 <- tryCatch({
  REEMctree_function(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = ~1 | Subject,
    correlation = nlme::corCAR1(form = ~1 | Subject), 
    ctree.control = party::ctree_control(mincriterion = 0.85),
    method = "REML" 
  )
}, error = function(e) {
  print("Error running REEMctree (Without Nested - Continuous AR(1)):")
  print(e)
  NULL
})
system_time_ctree_wn_car1 <- Sys.time() - start_time_ctree_wn_car1

metrics_ctree_wn_car1 <- list(MSE = NA, RMSE = NA, MAE = NA)
if (!is.null(reemctree_model_wn_car1) && !is.null(reemctree_model_wn_car1$EffectModel)) {
  residuals_wn_car1 <- reemctree_model_wn_car1$residuals[reemctree_model_wn_car1$Subset]
  metrics_ctree_wn_car1$MSE <- mean(residuals_wn_car1^2, na.rm = TRUE)
  metrics_ctree_wn_car1$RMSE <- sqrt(metrics_ctree_wn_car1$MSE)
  metrics_ctree_wn_car1$MAE <- mean(abs(residuals_wn_car1), na.rm = TRUE)
}
print(metrics_ctree_wn_car1)
print(system_time_ctree_wn_car1)


print("--- Starting REEMctree model: With Nested - AR(1) ---")
start_time_ctree_n_ar1 <- Sys.time()
reemctree_model_n_ar1 <- tryCatch({
  REEMctree_function(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = ~1 | Subject_Indices,
    correlation = nlme::corAR1(form = ~1 | Subject_Indices), 
    ctree.control = party::ctree_control(mincriterion = 0.85),
    method = "REML"
  )
}, error = function(e) {
  print("Error running REEMctree (With Nested - AR(1)):")
  print(e)
  NULL
})
system_time_ctree_n_ar1 <- Sys.time() - start_time_ctree_n_ar1

metrics_ctree_n_ar1 <- list(MSE = NA, RMSE = NA, MAE = NA)
if (!is.null(reemctree_model_n_ar1) && !is.null(reemctree_model_n_ar1$EffectModel)) {
  residuals_n_ar1 <- reemctree_model_n_ar1$residuals[reemctree_model_n_ar1$Subset]
  metrics_ctree_n_ar1$MSE <- mean(residuals_n_ar1^2, na.rm = TRUE)
  metrics_ctree_n_ar1$RMSE <- sqrt(metrics_ctree_n_ar1$MSE)
  metrics_ctree_n_ar1$MAE <- mean(abs(residuals_n_ar1), na.rm = TRUE)
}
print(metrics_ctree_n_ar1)
print(system_time_ctree_n_ar1)



table6_your_results <- data.frame(
  Method_Structure = character(),
  Correlation_Type = character(),
  RMSE = numeric(),
  MSE = numeric(),
  MAE = numeric(),
  System_Time_sec = numeric(),
  stringsAsFactors = FALSE
)


if(exists("metrics_ctree_simple") && !is.null(metrics_ctree_simple$RMSE) && !is.na(metrics_ctree_simple$RMSE)) {
  table6_your_results <- rbind(table6_your_results, data.frame(
    Method_Structure = "Unbiased RE-EM tree - Without nested structure",
    Correlation_Type = "Default",
    RMSE = metrics_ctree_simple$RMSE,
    MSE = metrics_ctree_simple$MSE,
    MAE = metrics_ctree_simple$MAE,
    System_Time_sec = as.numeric(system_time_ctree_simple, units = "secs")
  ))
}


if(!is.null(metrics_ctree_wn_car1$RMSE) && !is.na(metrics_ctree_wn_car1$RMSE)) {
  table6_your_results <- rbind(table6_your_results, data.frame(
    Method_Structure = "Unbiased RE-EM tree - Without nested structure",
    Correlation_Type = "Continuous AR(1)", 
    RMSE = metrics_ctree_wn_car1$RMSE,
    MSE = metrics_ctree_wn_car1$MSE,
    MAE = metrics_ctree_wn_car1$MAE,
    System_Time_sec = as.numeric(system_time_ctree_wn_car1, units = "secs")
  ))
}

if(exists("metrics_ctree_nested") && !is.null(metrics_ctree_nested$RMSE) && !is.na(metrics_ctree_nested$RMSE)) {
  table6_your_results <- rbind(table6_your_results, data.frame(
    Method_Structure = "Unbiased RE-EM tree - With nested structure",
    Correlation_Type = "Default",
    RMSE = metrics_ctree_nested$RMSE,
    MSE = metrics_ctree_nested$MSE,
    MAE = metrics_ctree_nested$MAE,
    System_Time_sec = as.numeric(system_time_ctree_nested, units = "secs")
  ))
}

if(!is.null(metrics_ctree_n_ar1$RMSE) && !is.na(metrics_ctree_n_ar1$RMSE)) {
  table6_your_results <- rbind(table6_your_results, data.frame(
    Method_Structure = "Unbiased RE-EM tree - With nested structure",
    Correlation_Type = "AR(1)",
    RMSE = metrics_ctree_n_ar1$RMSE,
    MSE = metrics_ctree_n_ar1$MSE,
    MAE = metrics_ctree_n_ar1$MAE,
    System_Time_sec = as.numeric(system_time_ctree_n_ar1, units = "secs")
  ))
}


if (nrow(table6_your_results) > 0) {
  table6_your_results_formatted <- table6_your_results
  table6_your_results_formatted$RMSE <- sprintf("%.1e", table6_your_results_formatted$RMSE)
  table6_your_results_formatted$MSE <- sprintf("%.1e", table6_your_results_formatted$MSE)
  table6_your_results_formatted$MAE <- sprintf("%.1e", table6_your_results_formatted$MAE)
  table6_your_results_formatted$System_Time_sec <- sprintf("%.2f", table6_your_results_formatted$System_Time_sec)
  
  cat("\n\n--- Performance metrics for Unbiased RE-EM tree (based on YOUR LogLik results) ---\n")
  print(table6_your_results_formatted, row.names = FALSE)
} else {
  cat("\n\n--- No valid model results to populate the performance metrics table. ---\n")
}
