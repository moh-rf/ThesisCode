
library(readr)
library(dplyr)
library(REEMtree)
library(nlme)
library(rpart)
library(ggplot2)
library(gridExtra) 
library(party)
library(partykit)

all_data <- read_csv("D:/payanname/Matlabcodes/HbO_Behavior_fnirs.csv")

data_filtered_tbl <- all_data %>%
  filter(n_back %in% c("2-back", "3-back"))

data_filtered_tbl$Subject <- as.factor(data_filtered_tbl$Subject)
data_filtered_tbl$Session <- as.factor(data_filtered_tbl$Session)
data_filtered_tbl$Gender <- as.factor(data_filtered_tbl$Gender)
data_filtered_tbl$n_back <- as.factor(data_filtered_tbl$n_back)
data_filtered_tbl$Indices <- as.factor(data_filtered_tbl$Indices)

data_filtered_df <- as.data.frame(data_filtered_tbl)


fixed_effects_formula <- value ~ MeanRT + n_back
random_effects_formula <- ~ 1 | Subject

# Model 1: RE-EM tree (simple random effects)

reem_model_simple_effects <- NULL 
system_time_reem_simple_effects <- NA
metrics_reem_simple_effects <- list(MSE = NA, RMSE = NA, MAE = NA)

start_time_simple_effects <- Sys.time() 
tryCatch({
  reem_model_simple_effects <- REEMtree(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = random_effects_formula)
  
  
  system_time_reem_simple_effects <- Sys.time() - start_time_simple_effects
  print("RE-EM tree model (simple random effects) built successfully.")
  
  predictions_simple_effects <- predict(reem_model_simple_effects, newdata = data_filtered_df) 
  actuals_simple_effects <- data_filtered_df$value 
  
  if (length(predictions_simple_effects) == length(actuals_simple_effects)) {
    residuals_simple_effects <- actuals_simple_effects - predictions_simple_effects 
    metrics_reem_simple_effects$MSE <- mean(residuals_simple_effects^2, na.rm = TRUE)
    metrics_reem_simple_effects$RMSE <- sqrt(metrics_reem_simple_effects$MSE)
    metrics_reem_simple_effects$MAE <- mean(abs(residuals_simple_effects), na.rm = TRUE)
  } else {
    print("Error: Length of predictions and actual values for the simple random effects model do not match.")
  }
  
}, error = function(e) {
  system_time_reem_simple_effects <- Sys.time() - start_time_simple_effects
  print("Error in building RE-EM tree model (simple random effects):")
  print(e)
})

print(paste("MSE (simple_effects):", format(metrics_reem_simple_effects$MSE, scientific = TRUE, digits = 2)))
print(paste("RMSE (simple_effects):", format(metrics_reem_simple_effects$RMSE, scientific = TRUE, digits = 2)))
print(paste("MAE (simple_effects):", format(metrics_reem_simple_effects$MAE, scientific = TRUE, digits = 2)))
print(paste("System Time (simple_effects, secs):", format(as.numeric(system_time_reem_simple_effects, units = "secs"), digits = 2)))
print("--- End of Model 1 ---")

#####################################################


data_filtered_df$Subject_Indices <- interaction(data_filtered_df$Subject, data_filtered_df$Indices)

random_effects_combined <- ~1 | Subject_Indices

# Start of Model 2
print("--- Start of Model 2: RE-EM tree (nested random effects using combined ID) ---")

reem_model_nested_effects <- NULL
system_time_reem_nested_effects <- NA
metrics_reem_nested_effects <- list(MSE = NA, RMSE = NA, MAE = NA)

start_time_nested_effects <- Sys.time()
tryCatch({
  reem_model_nested_effects <- REEMtree(
    formula = fixed_effects_formula,
    data = data_filtered_df,
    random = random_effects_combined
  )
  
  system_time_reem_nested_effects <- Sys.time() - start_time_nested_effects
  print("RE-EM tree model (nested random effects) built successfully.")
  
  predictions_nested_effects <- predict(reem_model_nested_effects, newdata = data_filtered_df)
  actuals_nested_effects <- data_filtered_df$value
  
  if (length(predictions_nested_effects) == length(actuals_nested_effects)) {
    residuals_nested_effects <- actuals_nested_effects - predictions_nested_effects
    metrics_reem_nested_effects$MSE <- mean(residuals_nested_effects^2, na.rm = TRUE)
    metrics_reem_nested_effects$RMSE <- sqrt(metrics_reem_nested_effects$MSE)
    metrics_reem_nested_effects$MAE <- mean(abs(residuals_nested_effects), na.rm = TRUE)
  } else {
    print("Error: Length of predictions and actual values for the nested model do not match.")
  }
  
}, error = function(e) {
  system_time_reem_nested_effects <- Sys.time() - start_time_nested_effects
  print("Error in building RE-EM tree model (nested random effects):")
  print(e)
})

# Display outputs
print(paste("MSE (nested_effects):", format(metrics_reem_nested_effects$MSE, scientific = TRUE, digits = 2)))
print(paste("RMSE (nested_effects):", format(metrics_reem_nested_effects$RMSE, scientific = TRUE, digits = 2)))
print(paste("MAE (nested_effects):", format(metrics_reem_nested_effects$MAE, scientific = TRUE, digits = 2)))
print(paste("System Time (nested_effects, secs):", format(as.numeric(system_time_reem_nested_effects, units = "secs"), digits = 2)))
print("--- End of Model 2 ---")


#####################################################

# --- Create dataframe for plotting Model 1 (simple) ---
plot_data_reem_simple <- NULL 
if (!is.null(reem_model_simple_effects) && exists("predictions_simple_effects") && exists("actuals_simple_effects")) {
  if (length(predictions_simple_effects) == length(actuals_simple_effects)) {
    valid_indices_simple <- !is.na(predictions_simple_effects) & !is.na(actuals_simple_effects)
    if (sum(valid_indices_simple) > 0) {
      plot_data_reem_simple <- data.frame(
        Fitted = predictions_simple_effects[valid_indices_simple],
        Actual = actuals_simple_effects[valid_indices_simple]
      )
    }
  }
}

# --- Create dataframe for plotting Model 2 (nested) ---
plot_data_reem_nested <- NULL
if (!is.null(reem_model_nested_effects) && exists("predictions_nested_effects") && exists("actuals_nested_effects")) {
  if (length(predictions_nested_effects) == length(actuals_nested_effects)) {
    valid_indices_nested <- !is.na(predictions_nested_effects) & !is.na(actuals_nested_effects)
    if (sum(valid_indices_nested) > 0) {
      plot_data_reem_nested <- data.frame(
        Fitted = predictions_nested_effects[valid_indices_nested],
        Actual = actuals_nested_effects[valid_indices_nested]
      )
    }
  }
}


plot_reem_a_final <- ggplot() + theme_void() + labs(title = "RE-EM Nested: Data Missing")
if (!is.null(plot_data_reem_nested) && nrow(plot_data_reem_nested) > 0) {
  plot_reem_a_final <- ggplot(plot_data_reem_nested, aes(x = Fitted, y = Actual)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(x = "Fitted Values", y = "Actual Values", 
         caption = "(a) RE-EM tree with nested data structure") + 
    theme_bw() + 
    theme(
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 10, b = 5)), 
      plot.margin = margin(t = 5, r = 5, b = 20, l = 5), 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
      axis.line = element_line(colour = "black", linewidth = 0.5) 
    ) +
    coord_cartesian(xlim = c(-0.002, 0.0025), ylim = c(-0.005, 0.0065), expand = TRUE) + 
    scale_x_continuous(breaks = seq(-0.002, 0.002, by = 0.001)) +
    scale_y_continuous(breaks = seq(-0.005, 0.005, by = 0.005))
}

plot_reem_b_final <- ggplot() + theme_void() + labs(title = "RE-EM Simple: Data Missing")
if (!is.null(plot_data_reem_simple) && nrow(plot_data_reem_simple) > 0) {
  plot_reem_b_final <- ggplot(plot_data_reem_simple, aes(x = Fitted, y = Actual)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(x = "Fitted Values", y = "Actual Values",
         caption = "(b) RE-EM tree without nested data structure") +
    theme_bw() + 
    theme(
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 10, b = 5)), 
      plot.margin = margin(t = 5, r = 5, b = 20, l = 5), 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
      axis.line = element_line(colour = "black", linewidth = 0.5) 
    ) +
    coord_cartesian(xlim = c(-0.0015, 0.0015), ylim = c(-0.005, 0.0065), expand = TRUE) +
    scale_x_continuous(breaks = seq(-0.001, 0.001, by = 0.001)) +
    scale_y_continuous(breaks = seq(-0.005, 0.005, by = 0.005))
}

grid.arrange(
    plot_reem_a_final, 
    plot_reem_b_final, 
    ncol = 2)

######################################## without nested structure

cor_structures <- list(
  "Default" = NULL,
  "AR1" = corAR1(form = ~ 1 | Subject),
  "CAR1" = corCAR1(form = ~ 1 | Subject),
  "CompSymm" = corCompSymm(form = ~ 1 | Subject)
)

aic_results_simple <- data.frame(Correlation = character(), AIC = numeric(), stringsAsFactors = FALSE)

for (name in names(cor_structures)) {
  cat("Fitting model with:", name, "\n")
  model <- tryCatch({
    REEMtree(
      formula = value ~ MeanRT + n_back,
      data = data_filtered_df,
      random = ~ 1 | Subject,
      correlation = cor_structures[[name]]
    )
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    aic_val <- AIC(model)
    aic_results_simple <- rbind(aic_results_simple, data.frame(Correlation = name, AIC = aic_val))
  } else {
    aic_results_simple <- rbind(aic_results_simple, data.frame(Correlation = name, AIC = NA))
  }
}

print(aic_results_simple)

###################################### with nested structure
cor_structures_nested <- list(
  "Default" = NULL,
  "AR1" = corAR1(form = ~ 1 | Subject_Indices),
  "CAR1" = corCAR1(form = ~ 1 | Subject_Indices),
  "CompSymm" = corCompSymm(form = ~ 1 | Subject_Indices)
)

aic_results_nested <- data.frame(Correlation = character(), AIC = numeric(), stringsAsFactors = FALSE)

for (name in names(cor_structures_nested)) {
  cat("Fitting nested model with:", name, "\n")
  model <- tryCatch({
    REEMtree(
      formula = value ~ MeanRT + n_back,
      data = data_filtered_df,
      random = ~ 1 | Subject_Indices,
      correlation = cor_structures_nested[[name]]
    )
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    aic_val <- AIC(model)
    aic_results_nested <- rbind(aic_results_nested, data.frame(Correlation = name, AIC = aic_val))
  } else {
    aic_results_nested <- rbind(aic_results_nested, data.frame(Correlation = name, AIC = NA))
  }
}

print(aic_results_nested)

#########################################################


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
  formula = value ~ MeanRT + n_back,
  data = data_filtered_df,
  random = ~1 | Subject
)
m1 <- compute_metrics(model1, data_filtered_df, data_filtered_df$value)
results_table <- rbind(results_table, data.frame(Structure = "Without nested", Correlation = "Default", t(m1)))

# Without nested - CAR(1) 
model2 <- REEMtree(
  formula = value ~ MeanRT + n_back,
  data = data_filtered_df,
  random = ~1 | Subject,
  correlation = corCAR1(form = ~1 | Subject)
)
m2 <- compute_metrics(model2, data_filtered_df, data_filtered_df$value)
results_table <- rbind(results_table, data.frame(Structure = "Without nested", Correlation = "CAR(1)", t(m2)))

# -----------------------
# With nested - Default
data_filtered_df$Subject_Indices <- interaction(data_filtered_df$Subject, data_filtered_df$Indices)

model3 <- REEMtree(
  formula = value ~ MeanRT + n_back,
  data = data_filtered_df,
  random = ~1 | Subject_Indices
)
m3 <- compute_metrics(model3, data_filtered_df, data_filtered_df$value)
results_table <- rbind(results_table, data.frame(Structure = "With nested", Correlation = "Default", t(m3)))

# With nested - AR(1)
model4 <- REEMtree(
  formula = value ~ MeanRT + n_back,
  data = data_filtered_df,
  random = ~1 | Subject_Indices,
  correlation = corAR1(form = ~1 | Subject_Indices)
)
m4 <- compute_metrics(model4, data_filtered_df, data_filtered_df$value)
results_table <- rbind(results_table, data.frame(Structure = "With nested", Correlation = "AR(1)", t(m4)))


print_scientific_table <- function(df, digits = 2) {
  df_formatted <- df
  df_formatted$RMSE <- formatC(df$RMSE, format = "e", digits = digits)
  df_formatted$MSE  <- formatC(df$MSE,  format = "e", digits = digits)
  df_formatted$MAE  <- formatC(df$MAE,  format = "e", digits = digits)
  print(df_formatted, row.names = FALSE)
}


print_scientific_table(results_table, digits = 2)

############################################


