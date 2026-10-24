############################################################
# Assignment 2 - Reproduction of missForest
# Paper: Stekhoven & Buhlmann - missForest
# Experiment: Parkinson continuous dataset
############################################################


# ==========================================================
# 1. LOAD PACKAGES
# ==========================================================

library(missForest)
library(imputeR)


# ==========================================================
# 2. LOAD DATA
# ==========================================================

data(parkinson, package = "imputeR")

# Store the original complete dataset
df <- parkinson

# Basic checks
print(dim(df))
print(sum(is.na(df)))

# Stop if the original dataset already contains missing values
stopifnot(sum(is.na(df)) == 0)


# ==========================================================
# 3. EXPERIMENT SETTINGS
# ==========================================================

# Missing-data percentages used in the experiment
missing_percentages <- c(0.10, 0.20, 0.30)

# Number of independent simulations per missing-data percentage
n_runs <- 50

# Number of trees in each Random Forest
ntree_value <- 100

# Base seed for reproducibility
base_seed <- 123


# ==========================================================
# 4. FUNCTION FOR ONE SIMULATION
# ==========================================================

run_simulation <- function(data_complete,
                           missing_percentage,
                           run_number,
                           ntree = 100,
                           base_seed = 123) {
  
  # --------------------------------------------------------
  # A. Set a reproducible but different seed for every run
  # --------------------------------------------------------
  
  current_seed <- base_seed +
    round(missing_percentage * 1000) +
    run_number
  
  set.seed(current_seed)
  
  
  # --------------------------------------------------------
  # B. Make a copy of the complete dataset
  # --------------------------------------------------------
  
  data_missing <- data_complete
  
  
  # --------------------------------------------------------
  # C. Determine how many cells must become missing
  # --------------------------------------------------------
  
  n_total <- nrow(data_missing) * ncol(data_missing)
  
  n_missing <- round(
    missing_percentage * n_total
  )
  
  
  # --------------------------------------------------------
  # D. Randomly select cells
  # --------------------------------------------------------
  
  missing_indices <- sample(
    seq_len(n_total),
    size = n_missing,
    replace = FALSE
  )
  
  missing_positions <- arrayInd(
    missing_indices,
    .dim = dim(data_missing)
  )
  
  
  # --------------------------------------------------------
  # E. Store the TRUE values before removing them
  # --------------------------------------------------------
  
  true_values <- data_complete[missing_positions]
  
  
  # --------------------------------------------------------
  # F. Replace selected values with NA
  # --------------------------------------------------------
  
  data_missing[missing_positions] <- NA
  
  
  # Safety check
  stopifnot(
    sum(is.na(data_missing)) == n_missing
  )
  
  
  # --------------------------------------------------------
  # G. Run missForest
  # --------------------------------------------------------
  
  start_time <- Sys.time()
  
  mf_result <- missForest(
    data_missing,
    ntree = ntree,
    verbose = FALSE
  )
  
  end_time <- Sys.time()
  
  
  # --------------------------------------------------------
  # H. Check that missForest filled all missing values
  # --------------------------------------------------------
  
  stopifnot(
    sum(is.na(mf_result$ximp)) == 0
  )
  
  
  # --------------------------------------------------------
  # I. Retrieve the imputed values
  # --------------------------------------------------------
  
  imputed_values <- mf_result$ximp[missing_positions]
  
  
  # --------------------------------------------------------
  # J. Calculate TRUE NRMSE
  #
  # This compares the imputed values with the original
  # values that we deliberately removed.
  # --------------------------------------------------------
  
  nrmse <- sqrt(
    mean((true_values - imputed_values)^2) /
      var(true_values)
  )
  
  
  # --------------------------------------------------------
  # K. Store OOB error reported by missForest
  # --------------------------------------------------------
  
  oob_error <- as.numeric(
    mf_result$OOBerror["NRMSE"]
  )
  
  
  # --------------------------------------------------------
  # L. Runtime
  # --------------------------------------------------------
  
  runtime_seconds <- as.numeric(
    difftime(
      end_time,
      start_time,
      units = "secs"
    )
  )
  
  
  # --------------------------------------------------------
  # M. Return results of this simulation
  # --------------------------------------------------------
  
  return(
    data.frame(
      missing_percentage = missing_percentage,
      run = run_number,
      seed = current_seed,
      n_missing = n_missing,
      NRMSE = nrmse,
      OOB_NRMSE = oob_error,
      runtime_seconds = runtime_seconds
    )
  )
}


# ==========================================================
# 5. RUN THE COMPLETE EXPERIMENT
# ==========================================================

all_results <- data.frame()

for (missing_percentage in missing_percentages) {
  
  cat(
    "\n============================================\n",
    "Missing percentage:",
    missing_percentage * 100,
    "%\n",
    "============================================\n"
  )
  
  for (run_number in 1:n_runs) {
    
    cat(
      "Running simulation",
      run_number,
      "of",
      n_runs,
      "...\n"
    )
    
    result <- run_simulation(
      data_complete = df,
      missing_percentage = missing_percentage,
      run_number = run_number,
      ntree = ntree_value,
      base_seed = base_seed
    )
    
    all_results <- rbind(
      all_results,
      result
    )
  }
}


# ==========================================================
# 6. INSPECT ALL INDIVIDUAL RESULTS
# ==========================================================

print(all_results)


# ==========================================================
# 7. CREATE SUMMARY TABLE
# ==========================================================

summary_results <- do.call(
  rbind,
  lapply(
    split(
      all_results,
      all_results$missing_percentage
    ),
    function(x) {
      
      data.frame(
        Missing_Percentage =
          unique(x$missing_percentage) * 100,
        
        Runs =
          nrow(x),
        
        Mean_NRMSE =
          mean(x$NRMSE),
        
        SD_NRMSE =
          sd(x$NRMSE),
        
        Mean_OOB_NRMSE =
          mean(x$OOB_NRMSE),
        
        SD_OOB_NRMSE =
          sd(x$OOB_NRMSE),
        
        Mean_Runtime_Seconds =
          mean(x$runtime_seconds)
      )
    }
  )
)

rownames(summary_results) <- NULL


# ==========================================================
# 8. PRINT FINAL RESULTS
# ==========================================================

cat(
  "\n\n============================================\n",
  "FINAL REPRODUCTION RESULTS\n",
  "============================================\n"
)

print(summary_results)


# ==========================================================
# 9. ADD MINIMUM AND MAXIMUM NRMSE
# ==========================================================

range_results <- do.call(
  rbind,
  lapply(
    split(
      all_results,
      all_results$missing_percentage
    ),
    function(x) {
      
      data.frame(
        Missing_Percentage =
          unique(x$missing_percentage) * 100,
        
        Minimum_NRMSE =
          min(x$NRMSE),
        
        Maximum_NRMSE =
          max(x$NRMSE)
      )
    }
  )
)

rownames(range_results) <- NULL

print(range_results)


# ==========================================================
# 10. CHECK OOB ERROR AGAINST TRUE ERROR
# ==========================================================

all_results$OOB_difference <- (
  all_results$OOB_NRMSE -
    all_results$NRMSE
)

all_results$OOB_absolute_difference <- abs(
  all_results$OOB_difference
)

oob_summary <- do.call(
  rbind,
  lapply(
    split(
      all_results,
      all_results$missing_percentage
    ),
    function(x) {
      
      data.frame(
        Missing_Percentage =
          unique(x$missing_percentage) * 100,
        
        Mean_True_NRMSE =
          mean(x$NRMSE),
        
        Mean_OOB_NRMSE =
          mean(x$OOB_NRMSE),
        
        Mean_Absolute_Difference =
          mean(x$OOB_absolute_difference)
      )
    }
  )
)

rownames(oob_summary) <- NULL

print(oob_summary)


# ==========================================================
# 11. SAVE RESULTS
# ==========================================================

write.csv(
  all_results,
  "missForest_all_runs.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "missForest_summary_results.csv",
  row.names = FALSE
)

write.csv(
  oob_summary,
  "missForest_OOB_comparison.csv",
  row.names = FALSE
)


# ==========================================================
# 12. SAVE SESSION INFORMATION
#
# Important for reproducibility:
# records R version and package versions.
# ==========================================================

sink("session_info.txt")
sessionInfo()
sink()


# ==========================================================
# 13. FINAL PARAMETER OVERVIEW
# ==========================================================

cat(
  "\nPARAMETERS USED\n",
  "Dataset: Parkinson\n",
  "Observations:", nrow(df), "\n",
  "Variables:", ncol(df), "\n",
  "Missing percentages: 10%, 20%, 30%\n",
  "Runs per percentage:", n_runs, "\n",
  "Number of trees:", ntree_value, "\n",
  "Base seed:", base_seed, "\n",
  "Evaluation metric: NRMSE\n"
)


getwd()
list.files()
