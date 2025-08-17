
library(gpboost)
library(dplyr)
library(readr)
library(knitr)
library(ggplot2)
library(patchwork)


set.seed(123) 

# STEP 1: DATA PREPARATION 
file_path <- "D:\\payanname\\Matlabcodes\\HbO_Behavior_fnirs.csv"
all_data <- read_csv(file_path)

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
y_original <- data$value


SCALING_FACTOR <- 1e4
y_scaled <- y_original * SCALING_FACTOR
dtrain <- gpb.Dataset(data = X, label = y_scaled)

results_table4 <- data.frame()

# ===================================================================
# MODEL 1: GPBoost without nested structure (using GD)
# ===================================================================
cat("Running GPBoost WITHOUT nested structure using GD optimizer...\n")

group_data_wn <- data$Subject
gp_model_wn <- GPModel(group_data = group_data_wn, likelihood = "gaussian")
gp_model_wn$set_optim_params(params = list(optimizer_cov = "gradient_descent"))

start_time <- proc.time()
bst_wn <- gpboost(
  data = dtrain, gp_model = gp_model_wn,
  nrounds = 150, learning_rate = 0.1, verbose = 0 
)
end_time <- proc.time()
exec_time <- (end_time - start_time)["user.self"]

preds_scaled <- predict(bst_wn, data = X, group_data_pred = group_data_wn)$response_mean
preds_unscaled <- preds_scaled / SCALING_FACTOR

mse <- mean((y_original - preds_unscaled)^2)
rmse <- sqrt(mse)
mae <- mean(abs(y_original - preds_unscaled))

results_table4 <- rbind(results_table4, data.frame(
  Method = "GPBoost without nested structure",
  RMSE = rmse, MSE = mse, MAE = mae, System_Time = exec_time
))

# ===================================================================
# MODEL 2: GPBoost with nested structure (using GD)
# ===================================================================
cat("Running GPBoost WITH nested structure using GD optimizer...\n")

group_data_n <- as.matrix(data[, c("Subject", "Indices")])
gp_model_n <- GPModel(group_data = group_data_n, likelihood = "gaussian")
gp_model_n$set_optim_params(params = list(optimizer_cov = "gradient_descent"))

start_time <- proc.time()
bst_n <- gpboost(
  data = dtrain, gp_model = gp_model_n,
  nrounds = 127, learning_rate = 0.1, verbose = 0 
)
end_time <- proc.time()
exec_time <- (end_time - start_time)["user.self"]

preds_scaled_n <- predict(bst_n, data = X, group_data_pred = group_data_n)$response_mean
preds_unscaled_n <- preds_scaled_n / SCALING_FACTOR

mse_n <- mean((y_original - preds_unscaled_n)^2)
rmse_n <- sqrt(mse_n)
mae_n <- mean(abs(y_original - preds_unscaled_n))

results_table4 <- rbind(results_table4, data.frame(
  Method = "GPBoost with nested structure",
  RMSE = rmse_n, MSE = mse_n, MAE = mae_n, System_Time = exec_time
))

# STEP 4: DISPLAY THE FINAL TABLE
formatted_results <- results_table4 %>%
  mutate(
    RMSE = format(RMSE, scientific = TRUE, digits = 3),
    MSE = format(MSE, scientific = TRUE, digits = 3),
    MAE = format(MAE, scientific = TRUE, digits = 3),
    System_Time = round(System_Time, 2)
  )

suppressWarnings({
  print(kable(formatted_results, format = "pipe", caption = "Final and Correct GPBoost Rows for Table 4 (using GD optimizer)"))
})

##################################################

# STEP 5: REPLICATE FIGURE 9 FROM THE APPENDIX

y <- data$value
# --- Create data frames for plotting ---
# The plots need data frames that contain both the actual and fitted values.

# Data for the "with nested structure" plot
plot_df_nested <- data.frame(
  Actual = y, # Original, unscaled response variable
  Fitted = preds_unscaled_n # Predictions from the nested model, rescaled
)

# Data for the "without nested structure" plot
plot_df_without <- data.frame(
  Actual = y, # Original, unscaled response variable
  Fitted = preds_unscaled # Predictions from the simple model, rescaled
)

# --- Create Plot (a): GPBoost with nested structure ---
plot_a <- ggplot(data = plot_df_nested, aes(x = Fitted, y = Actual)) +
  # Add the scatter plot points
  geom_point(alpha = 0.8) + 
  # Add a red line of best fit (linear model) to show the trend
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  # Set titles and labels
  labs(
    title = "Fitted vs. actual values for GPBoost with nested structure",
    x = "Fitted Values",
    y = "Actual Values",
    caption = "(a) GPBoost with nested data structure" # Add the bottom label
  ) +
  # Use a clean, professional theme
  theme_bw() +
  # Center the caption at the bottom
  theme(plot.caption = element_text(hjust = 0.5, size = 12))


# --- Create Plot (b): GPBoost without nested structure ---
plot_b <- ggplot(data = plot_df_without, aes(x = Fitted, y = Actual)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  labs(
    title = "Fitted vs. actual values for GPBoost without nested structure",
    x = "Fitted Values",
    y = "Actual Values",
    caption = "(b) GPBoost without nested data structure" # Add the bottom label
  ) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 12))


# --- Combine the two plots side-by-side using patchwork ---
final_figure <- plot_a + plot_b

# Display the final figure
print(final_figure)

####################################################

X <- model.matrix(value ~ Age + Accuracy + MeanRT + n_back + Session + Gender - 1, data = data)
y_original <- data$value
N <- length(y_original)

SCALING_FACTOR <- 1e4
y_scaled_label <- y_original * SCALING_FACTOR

dtrain <- gpb.Dataset(data = X, label = y_scaled_label)

LEARNING_RATE_ARTICLE <- 0.1

# =================================================================================
cat("===== Running Models WITHOUT Nested Structure (Article Simulation) =====\n")

results_without_nested <- data.frame()
group_data_wn <- data$Subject
optimizers_wn <- list(
  GD = list(name = "gradient_descent", nrounds = 150),
  WLS = list(name = "newton", nrounds = 4),
  NM = list(name = "nelder_mead", nrounds = 2),
  BFGS = list(name = "lbfgs", nrounds = 3)
)

for (opt_name in names(optimizers_wn)) {
  cat(paste("Running Optimizer:", opt_name, "...\n"))
  
  params <- optimizers_wn[[opt_name]]
  gp_model <- GPModel(group_data = group_data_wn, likelihood = "gaussian")
  gp_model$set_optim_params(params = list(optimizer_cov = params$name))
  
  start_time <- proc.time()
  bst <- gpboost(
    data = dtrain, gp_model = gp_model, nrounds = params$nrounds,
    learning_rate = LEARNING_RATE_ARTICLE, verbose = 0, valids = list(train = dtrain), record = TRUE
  )
  end_time <- proc.time()
  exec_time <- (end_time - start_time)["user.self"]
  
  preds_scaled <- predict(bst, data = X, group_data_pred = group_data_wn)$response_mean
  preds_unscaled <- preds_scaled / SCALING_FACTOR
  
  mse <- mean((y_original - preds_unscaled)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_original - preds_unscaled))
  
  nll_scaled <- tail(gpb.get.eval.result(bst, "train", "Negative log-likelihood"), 1)
  nll_original <- nll_scaled - N * log(SCALING_FACTOR)
  
  results_without_nested <- rbind(results_without_nested, data.frame(
    Optimizer = opt_name, RMSE = rmse, MSE = mse, MAE = mae,
    Iterations = params$nrounds, System_Time = round(exec_time, 2), Neg_Log_Likelihood = nll_original
  ))
}

# =================================================================================
cat("\n===== Running Models WITH Nested Structure (Article Simulation) =====\n")

results_with_nested <- data.frame()
group_data_n <- as.matrix(data[, c("Subject", "Indices")])
optimizers_n <- list(
  GD = list(name = "gradient_descent", nrounds = 127),
  WLS = list(name = "newton", nrounds = 6),
  NM = list(name = "nelder_mead", nrounds = 2),
  BFGS = list(name = "lbfgs", nrounds = 2)
)

for (opt_name in names(optimizers_n)) {
  cat(paste("Running Optimizer:", opt_name, "...\n"))
  
  params <- optimizers_n[[opt_name]]
  gp_model <- GPModel(group_data = group_data_n, likelihood = "gaussian")
  gp_model$set_optim_params(params = list(optimizer_cov = params$name))
  
  start_time <- proc.time()
  bst <- gpboost(
    data = dtrain, gp_model = gp_model, nrounds = params$nrounds,
    learning_rate = LEARNING_RATE_ARTICLE, verbose = 0, valids = list(train = dtrain), record = TRUE
  )
  end_time <- proc.time()
  exec_time <- (end_time - start_time)["user.self"]
  
  preds_scaled <- predict(bst, data = X, group_data_pred = group_data_n)$response_mean
  preds_unscaled <- preds_scaled / SCALING_FACTOR
  
  mse <- mean((y_original - preds_unscaled)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_original - preds_unscaled))
  
  nll_scaled <- tail(gpb.get.eval.result(bst, "train", "Negative log-likelihood"), 1)
  nll_original <- nll_scaled - N * log(SCALING_FACTOR)
  
  results_with_nested <- rbind(results_with_nested, data.frame(
    Optimizer = opt_name, RMSE = rmse, MSE = mse, MAE = mae,
    Iterations = params$nrounds, System_Time = round(exec_time, 2), Neg_Log_Likelihood = nll_original
  ))
}

# STEP 3
format_results_full <- function(df) {
  df %>%
    mutate(
      RMSE = format(RMSE, scientific = TRUE, digits = 3),
      MSE = format(MSE, scientific = TRUE, digits = 3),
      MAE = format(MAE, scientific = TRUE, digits = 3),
      Neg_Log_Likelihood = round(Neg_Log_Likelihood, 1)
    )
}

cat("\n\n--- Table 5 Recreation (Article Simulation): Without Nested Structure ---\n")
suppressWarnings(print(kable(format_results_full(results_without_nested), format = "pipe",
                             col.names = c("Optimizer", "RMSE", "MSE", "MAE", "Iterations", "System Time", "Negative log-likelihood"))))

cat("\n\n--- Table 5 Recreation (Article Simulation): With Nested Structure ---\n")
suppressWarnings(print(kable(format_results_full(results_with_nested), format = "pipe",
                             col.names = c("Optimizer", "RMSE", "MSE", "MAE", "Iterations", "System Time", "Negative log-likelihood"))))

