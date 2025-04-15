#!/usr/bin/env Rscript
# summarize_results.R
#
# This script takes in one or more RData files (each containing a data frame
# named final_results with simulation output) and produces a summary data frame.
# Each summary row is summarized for the simulation output from that run and
# includes the mean and standard deviation for the columns:
#   - mse_rInstance
#   - mse_Dapper
#   - avg_post_means
#   - avg_post_vars
#
# Usage:
#   Rscript summarize_results.R /path/to/run1/sim_results.RData /path/to/run2/sim_results.RData ...

# Get command line arguments (the file paths)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please supply one or more RData file paths as command line arguments.")
}

# Load necessary package
suppressPackageStartupMessages(library(dplyr))

# Initialize an empty list to hold summary rows
results_list <- list()

# Loop over each file passed as an argument
for (file in args) {
  
  # Check if the file exists
  if (!file.exists(file)) {
    warning(sprintf("File '%s' does not exist. Skipping...", file))
    next
  }
  
  # Load the RData file. It is assumed that it loads a data frame called final_results.
  load(file)  # final_results should now be in the workspace
  
  # Basic check that final_results exists and is a data frame
  if (!exists("final_results") || !is.data.frame(final_results)) {
    warning(sprintf("File '%s' did not load a valid final_results data frame. Skipping...", file))
    next
  }
  
  # Extract the run name (parameter combo) from the file path.
  # For example, if file is 
  # "/home/navarr72/dpwork/dapper_comparison_output_1000_100_4_1_4_2/sim_results.RData",
  # then run_name will be "dapper_comparison_output_1000_100_4_1_4_2".
  run_name <- basename(dirname(file))
  
  # Compute summary statistics over replicates in final_results.
  summary_stats <- final_results %>%
    summarise(
      mse_rInstance_mean = mean(mse_rInstance, na.rm = TRUE),
      mse_rInstance_sd   = sd(mse_rInstance, na.rm = TRUE),
      mse_Dapper_mean    = mean(mse_Dapper, na.rm = TRUE),
      mse_Dapper_sd      = sd(mse_Dapper, na.rm = TRUE),
      avg_post_means_mean = mean(avg_post_means, na.rm = TRUE),
      avg_post_means_sd   = sd(avg_post_means, na.rm = TRUE),
      avg_post_vars_mean  = mean(avg_post_vars, na.rm = TRUE),
      avg_post_vars_sd    = sd(avg_post_vars, na.rm = TRUE)
    ) %>%
    # Temporarily keep the run_name for setting row names later
    mutate(run_name = run_name)
  
  # Add the summary to our results list
  results_list[[run_name]] <- summary_stats
  
  # Remove final_results to avoid conflict on the next iteration
  rm(final_results)
}

# Combine all summaries into one data frame
final_summary_df <- bind_rows(results_list)

# Set row names using the run_name values, then remove the run_name column so it only appears as row names.
rownames(final_summary_df) <- final_summary_df$run_name
final_summary_df <- final_summary_df %>% select(-run_name)

# Print the summary data frame to the console
print(final_summary_df)

# Save the summary results to an RData file
output_rdata <- "summary_results.RData"
save(final_summary_df, file = output_rdata)
cat(sprintf("Summary results saved to %s\n", output_rdata))
