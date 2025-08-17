
# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(nlme)
library(gridExtra)

# Load the data
df <- read.csv("D:/payanname/Matlabcodes/HbO_Behavior_fnirs.csv")

# Convert variables to factors as described in the paper or implied
df$n_back <- factor(df$n_back, levels = c("0-back", "2-back", "3-back"), labels = c("0", "2", "3")) # Use simpler labels if needed for model formulas
df$Subject <- as.factor(df$Subject)
df$Indices <- as.factor(df$Indices)
df$Session <- as.factor(df$Session) 
df$Gender <- as.factor(df$Gender)  

# Filter for 2-back and 3-back conditions as per the paper
df_filtered <- df %>%
  filter(n_back %in% c("2", "3")) %>% # Use new labels if changed above
  droplevels() # Drop unused factor levels after filtering

# Define fixed and random effects based on the parsimonious model in the paper
# The paper states: "We have n-back and MeanRT variables in the parsimonious model."
fixed_effects_formula <- value ~ n_back + MeanRT
random_effects_formula <- ~ 1 | Subject/Indices

# Control parameters for lme, especially for potentially complex models
ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", msVerbose = FALSE)

# Model 1: No within-group correlation (Default)
# Argument in paper: Default
cat("Fitting Model 1 (No correlation)...\n")
m1 <- lme(fixed = fixed_effects_formula,
          random = random_effects_formula,
          data = df_filtered,
          method = "ML", # Use ML for AIC comparison
          control = ctrl)
aic_m1 <- AIC(m1)

# Model 2: Autoregressive process of order 1 (AR(1))
# Argument in paper: corAR1()
cat("Fitting Model 2 (AR(1))...\n")
m2 <- lme(fixed = fixed_effects_formula,
          random = random_effects_formula,
          correlation = corAR1(), 
          data = df_filtered,
          method = "ML",
          control = ctrl)
aic_m2 <- AIC(m2)

# Model 3: Continuous AR(1)
# Argument in paper: corCAR1()
cat("Fitting Model 3 (Continuous AR(1))...\n")
m3 <- lme(fixed = fixed_effects_formula,
          random = random_effects_formula,
          correlation = corCAR1(), 
          data = df_filtered,
          method = "ML",
          control = ctrl)
aic_m3 <- AIC(m3)

# Model 4: A general positive definite matrix (Unstructured/General Symmetric)
# Argument in paper: corSymm()
# This model can be very complex if Indices has many levels.
cat("Fitting Model 4 (General Symmetric - corSymm)...\n")
aic_m4 <- NA # Initialize
tryCatch({
  m4 <- lme(fixed = fixed_effects_formula,
            random = random_effects_formula,
            # corSymm() without form applies to the innermost grouping factor (Indices within Subject)
            # using its levels as "occasions".
            correlation = corSymm(form = ~ 1 | Subject/Indices), # Explicitly stating for clarity
            data = df_filtered,
            method = "ML",
            control = lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim", niterEM = 100, msVerbose = TRUE)) # Increased iterations
  aic_m4 <- AIC(m4)
  if (m4$numIter > ctrl$maxIter || !m4$apVar) warning("Model m4 (corSymm) may have convergence issues despite producing an AIC.")
}, error = function(e) {
  warning("Model m4 (corSymm) failed to converge or other error: ", e$message)
  # m4 might not be assigned or might be problematic, so AIC remains NA
})


# Model 5: A compound-symmetric matrix
# Argument in paper: corCompSymm()
cat("Fitting Model 5 (Compound Symmetry - corCompSymm)...\n")
m5 <- lme(fixed = fixed_effects_formula,
          random = random_effects_formula,
          # corCompSymm() without form applies to the innermost grouping factor
          correlation = corCompSymm(form = ~ 1 | Subject/Indices), # Explicitly stating for clarity
          data = df_filtered,
          method = "ML",
          control = ctrl)
aic_m5 <- AIC(m5)

# Create AIC table
aic_table <- data.frame(
  Model_Description = c("No within-group correlation (Default)",
                        "Autoregressive process of order 1 (AR(1))",
                        "Continuous AR(1)",
                        "A general positive definite matrix (corSymm)",
                        "A compound-symmetric matrix (corCompSymm)"),
  Argument_Paper = c("Default", "corAR1()", "corCAR1()", "corSymm()", "corCompSymm()"),
  AIC_Calculated = c(aic_m1, aic_m2, aic_m3, aic_m4, aic_m5)
)

print(aic_table)

########################################################

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




plot_data_no_nest <- data.frame(Fitted = predicted_lmm_no_nest, Actual = actual_values, Model = "LMM without nested structure")
plot_data_nest <- data.frame(Fitted = predicted_m4, Actual = actual_values, Model = "LMM with nested structure (corSymm)")

#LMM without nested structure 
p_b <- ggplot(plot_data_no_nest, aes(x = Fitted, y = Actual)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  labs(x = "Fitted Values", y = "Actual Values") +
  theme_bw()

#LMM with nested structure
p_a <- ggplot(plot_data_nest, aes(x = Fitted, y = Actual)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  labs(x = "Fitted Values", y = "Actual Values") +
  theme_bw()


label_a <- textGrob("(a) LMM with nested data structure", gp = gpar(fontsize = 10), hjust = 0.5)
label_b <- textGrob("(b) LMM without nested data structure", gp = gpar(fontsize = 10), hjust = 0.5)


grob_a <- arrangeGrob(p_a, label_a, ncol = 1, heights = c(10, 1))
grob_b <- arrangeGrob(p_b, label_b, ncol = 1, heights = c(10, 1))


cat("\nGenerating plots...\n")
grid.arrange(grob_a, grob_b, ncol = 2)



###########################################################################


# Fixed effects formula 
fixed_effects_formula <- value ~ n_back + MeanRT

# Random effects formula 
random_effects_nonnested <- ~ 1 | Subject


ctrl_simple <- lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", msVerbose = FALSE)


cat("Fitting Model 1 (No within-group correlation)...\n")
m1_nn <- lme(fixed = fixed_effects_formula,
             random = random_effects_nonnested,
             data = df_filtered,
             method = "ML",
             control = ctrl_simple)
aic_m1_nn <- AIC(m1_nn)

#AR(1)
cat("Fitting Model 2 (AR(1))...\n")
m2_nn <- lme(fixed = fixed_effects_formula,
             random = random_effects_nonnested,
             correlation = corAR1(form = ~1 | Subject),
             data = df_filtered,
             method = "ML",
             control = ctrl_simple)
aic_m2_nn <- AIC(m2_nn)

# Continuous AR(1)
cat("Fitting Model 3 (Continuous AR(1))...\n")
m3_nn <- lme(fixed = fixed_effects_formula,
             random = random_effects_nonnested,
             correlation = corCAR1(form = ~1 | Subject),
             data = df_filtered,
             method = "ML",
             control = ctrl_simple)
aic_m3_nn <- AIC(m3_nn)

#General positive definite (corSymm)
cat("Fitting Model 4 (General Symmetric - corSymm)...\n")
aic_m4_nn <- NA 
tryCatch({
  m4_nn <- lme(fixed = fixed_effects_formula,
               random = random_effects_nonnested,
               correlation = corSymm(form = ~1 | Subject),
               data = df_filtered,
               method = "ML",
               control = lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim", niterEM = 100))
  aic_m4_nn <- AIC(m4_nn)
}, error = function(e) {
  warning("Model m4_nn (corSymm) failed: ", e$message)
})

#Compound symmetry
cat("Fitting Model 5 (Compound Symmetry - corCompSymm)...\n")
m5_nn <- lme(fixed = fixed_effects_formula,
             random = random_effects_nonnested,
             correlation = corCompSymm(form = ~1 | Subject),
             data = df_filtered,
             method = "ML",
             control = ctrl_simple)
aic_m5_nn <- AIC(m5_nn)


aic_table_nonnested <- data.frame(
  Model_Description = c("No within-group correlation (Default)",
                        "Autoregressive process of order 1 (AR(1))",
                        "Continuous AR(1)",
                        "General positive definite matrix (corSymm)",
                        "Compound-symmetric matrix (corCompSymm)"),
  Structure = "Non-nested",
  AIC = c(aic_m1_nn, aic_m2_nn, aic_m3_nn, aic_m4_nn, aic_m5_nn)
)

cat("\n--- AIC Table (Non-nested models) ---\n")
print(aic_table_nonnested, row.names = FALSE)





