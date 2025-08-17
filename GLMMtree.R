
### LOAD LIBRARIES
library(glmertree)
library(dplyr)
library(readr)
library(lme4) 
library(ggplot2)
library(gridExtra)
library(partykit)

### STEP 1: DATA PREPARATION
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
cat("--- Step 2: Fitting Models and Measuring Time ---\n")

time_subj_only <- system.time({
  lt <- lmertree(value ~ n_back + MeanRT + Session + Accuracy | Subject | Age + Gender,
                 data = data_for_glmm)
})

glmm_tree_subj_only <- lt
print(lt)
plot(lt, which = "all")
plot(as.party(lt), tp_args = list(gp = gpar(fontsize = 1)))


time_nested <- system.time({
  glmm_tree_nested <- lmertree(value ~ n_back + MeanRT + Session + Accuracy | (1 | Subject/Indices) | Age + Gender,
                               data = data_for_glmm)
})
plot(glmm_tree_nested, which = "tree")


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


### STEP 4: DISPLAYING COMBINED RESULTS
cat("\n--- Step 4: Combined Results Table ---\n")


results_table <- data.frame(
  Method = c("GLMM tree without nested structure", "GLMM tree with nested structure"),
  RMSE = c(metrics_glmm_subj_only$RMSE, metrics_glmm_nested$RMSE),
  MSE = c(metrics_glmm_subj_only$MSE, metrics_glmm_nested$MSE),
  MAE = c(metrics_glmm_subj_only$MAE, metrics_glmm_nested$MAE),
  SystemTime_sec = c(time_subj_only["elapsed"], time_nested["elapsed"])
)


print(results_table, row.names = FALSE, digits = 4)


### STEP 4: VISUALIZATION – Fitted vs. Actual

# Predictions
predictions_subj_only <- predict(glmm_tree_subj_only, newdata = data_for_glmm, type = "response", allow.new.levels = TRUE)
predictions_nested <- predict(glmm_tree_nested, newdata = data_for_glmm, type = "response", allow.new.levels = TRUE)
actual_values_glmm <- data_for_glmm$value 

# Create plot dataframes
plot_data_glmm_subj_only <- na.omit(data.frame(Fitted = predictions_subj_only, Actual = actual_values_glmm))
plot_data_glmm_nested <- na.omit(data.frame(Fitted = predictions_nested, Actual = actual_values_glmm))

# Panel B – without nested structure
plot_glmm_b <- ggplot(plot_data_glmm_subj_only, aes(x = Fitted, y = Actual)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") + 
  labs(title = "GLMM tree without nested structure",
       x = "Fitted Values", y = "Actual Values",
       caption = "(b) GLMM tree without nested data structure" ) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 12))

# Panel A – with nested structure
plot_glmm_a <- ggplot(plot_data_glmm_nested, aes(x = Fitted, y = Actual)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") + 
  labs(title = "GLMM tree with nested structure",
       x = "Fitted Values", y = "Actual Values",
       caption = "(a) GLMM tree with nested data structure") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 12))

cat("\nGenerating Fitted vs. Actual plots...\n")

print(grid.arrange(plot_glmm_a, plot_glmm_b, ncol = 2, 
                   top = grid::textGrob("Fitted vs. Actual for GLMM Trees", 
                                        gp = grid::gpar(fontsize = 14))))

###########################################################
cat("\nGenerating plot similar to Figure 5: Fitted vs. Actual for GLMM tree (nested) for Subject-7 by Indices...\n")

selected_subject_id_fig5 <- "7" 
ncols_for_indices_fig5 <- 9    

data_subject_7_fig5 <- data_for_glmm %>%
  filter(Subject == selected_subject_id_fig5)


predictions_fig5 <- predict(glmm_tree_nested, newdata = data_subject_7_fig5, type = "response", allow.new.levels = TRUE)
actual_values_fig5 <- data_subject_7_fig5$value
indices_fig5 <- data_subject_7_fig5$Indices 

plot_df_fig5 <- data.frame(
  Fitted = predictions_fig5,
  Actual = actual_values_fig5,
  Indices = indices_fig5
)

plot_df_fig5 <- na.omit(plot_df_fig5)


plot_figure_5_replication <- ggplot(plot_df_fig5, aes(x = Fitted, y = Actual)) +
  geom_point(alpha = 0.7, size = 0.9, color = "black", shape = 16) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linewidth = 0.5, fullrange = TRUE) +
  facet_wrap(~ Indices, ncol = ncols_for_indices_fig5, scales = "fixed") +
  labs(
    title = paste("Figure 5. Fitted versus actual values plot of GLMM tree for indices nested within subject-", selected_subject_id_fig5, sep = ""),
    x = "Fitted Values",
    y = "Actual Values"
  ) +
  scale_x_continuous(
    breaks = c(-0.0015, -0.0005, 0.0005),
    labels = function(x) sprintf("%.4f", x)
    # limits = c(min_fitted_val - offset, max_fitted_val + offset)
  ) +
  scale_y_continuous(
    breaks = c(-0.003, -0.002, -0.001, 0.000, 0.001),
    labels = function(y) sprintf("%.3f", y)
    # limits = c(min_actual_val - offset, max_actual_val + offset)
  ) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "moccasin", color = "grey50", linewidth = 0.5),
    strip.text = element_text(size = 7, face = "bold", margin = margin(t = 2, b = 2)),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(margin = margin(t = 1)),
    panel.spacing = unit(0.2, "lines"),
    panel.grid.major = element_line(colour = "grey95"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
  )

cat("Plot similar to Figure 5 generated. Displaying now...\n")
print(plot_figure_5_replication)

############################################################################
