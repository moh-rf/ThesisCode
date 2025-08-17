
### LOAD LIBRARIES
library(LongCART)
library(dplyr)
library(readr)

### STEP 1: DATA PREPARATION
file_path <- "D:\\payanname\\Matlabcodes\\HbO_Behavior_fnirs.csv"
all_data <- read_csv(file_path, show_col_types = FALSE)

# Filter for 2-back and 3-back
data <- all_data %>%
  filter(n_back %in% c("2-back", "3-back")) %>%
  mutate(
    Subject = as.factor(Subject),
    Session = as.factor(Session),
    Gender = as.factor(Gender),
    n_back = factor(n_back, levels = c("2-back", "3-back")),
    Indices = as.factor(Indices)
  )

### STEP 2: DEFINE FIXED EFFECTS MODEL
fixed_formula <- as.formula("value ~ MeanRT + n_back")

### STEP 3: FIT LONGCART MODEL (WITH AGE AND GENDER)
longcart_model <- tryCatch({
  LongCART(
    data = data,
    patid = "Subject",
    fixed = fixed_formula,
    gvars = c("Gender", "Age"),    
    tgvars = c(0, 1),              
    minsplit = 40,
    minbucket = 20,
    alpha = 0.05,
    coef.digits = 2,
    print.lme = TRUE
  )
}, error = function(e) {
  cat("Error during LongCART execution:\n", e$message, "\n")
  NULL
})

# STEP 4: CHECK RESULTS
if (!is.null(longcart_model)) {
  summary(longcart_model)     # Print a summary of the fitted model
  plot(longcart_model)        # Plot the tree if any split has occurred
} else {
  cat("The model did not run or no splitting was performed.\n")
}

