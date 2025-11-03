## Analysis script for Conditional A DORM simulation results with L=5

cat("=== Initializing analysis ===\n")

# Setup paths
if (requireNamespace("here", quietly = TRUE)) {
  library(here)
} else {
  here <- function(...) file.path(getwd(), ...)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Define input/output directories
data_dir <- here('simu', 'simu_10242025', 'data')
output_dir <- here('simu', 'simu_10242025', 'results')

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

## Function to load and combine results across seeds
load_and_combine <- function(data_dir, pattern, result_type = "worst") {
  # Find all matching files
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  
  if(length(files) == 0) {
    stop("No files found matching pattern: ", pattern)
  }
  
  cat(sprintf("Found %d files for %s case\n", length(files), result_type))
  
  # Load all results
  all_results <- list()
  all_params <- list()
  
  for(i in seq_along(files)) {
    data <- readRDS(files[i])
    seed <- data$parameters$seed
    
    if(result_type == "worst") {
      df <- data$final_worst
    } else {
      df <- data$final_ave
    }
    
    df$seed <- seed
    all_results[[i]] <- df
    all_params[[i]] <- data$parameters
  }
  
  # Combine all results
  combined_df <- do.call(rbind, all_results)
  
  return(list(
    data = combined_df,
    params = all_params[[1]],  # Parameters should be same across seeds
    n_seeds = length(files)
  ))
}

## Load worst case results
cat("\n=== Loading worst case results ===\n")
worst_results <- load_and_combine(data_dir, "CondA_L5_worst_seed.*\\.rds$", "worst")

## Load average case results
cat("\n=== Loading average case results ===\n")
ave_results <- load_and_combine(data_dir, "CondA_L5_ave_seed.*\\.rds$", "average")

## Function to compute summary statistics across seeds
compute_summary <- function(df, group_vars = c("smax", "true-s")) {
  # Methods to summarize
  methods <- c("Ours", "SS", "SA", "RA", "MI", "PA")
  
  # Compute mean and standard error across seeds
  summary_df <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      across(all_of(methods), 
             list(mean = ~mean(.x, na.rm = TRUE),
                  se = ~sd(.x, na.rm = TRUE) / sqrt(n()),
                  sd = ~sd(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_seeds = n(),
      .groups = "drop"
    )
  
  return(summary_df)
}

## Compute summaries
cat("\n=== Computing summary statistics ===\n")
worst_summary <- compute_summary(worst_results$data)
ave_summary <- compute_summary(ave_results$data)

cat(sprintf("Summarized over %d seeds\n", worst_results$n_seeds))

## Save summary tables
write.csv(worst_summary, 
          file.path(output_dir, "worst_case_summary.csv"), 
          row.names = FALSE)
write.csv(ave_summary, 
          file.path(output_dir, "average_case_summary.csv"), 
          row.names = FALSE)

cat("Summary tables saved\n")

## Function to create performance plots
create_performance_plot <- function(summary_df, case_type = "Worst Case") {
  # Reshape data for plotting
  plot_data <- summary_df %>%
    select(smax, `true-s`, ends_with("_mean")) %>%
    pivot_longer(cols = ends_with("_mean"),
                 names_to = "method",
                 values_to = "error") %>%
    mutate(method = gsub("_mean", "", method))
  
  # Add standard errors for error bars
  se_data <- summary_df %>%
    select(smax, `true-s`, ends_with("_se")) %>%
    pivot_longer(cols = ends_with("_se"),
                 names_to = "method",
                 values_to = "se") %>%
    mutate(method = gsub("_se", "", method))
  
  plot_data <- left_join(plot_data, se_data, by = c("smax", "true-s", "method"))
  
  # Create plot for each smax value
  plots <- list()
  smax_values <- unique(plot_data$smax)
  
  for(i in seq_along(smax_values)) {
    smax_val <- smax_values[i]
    df_subset <- plot_data %>% filter(smax == smax_val)
    
    p <- ggplot(df_subset, aes(x = `true-s`, y = error, color = method, group = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = error - se, ymax = error + se), 
                    width = 0.02, alpha = 0.5) +
      labs(title = sprintf("%s (smax = %.2f)", case_type, smax_val),
           x = "True Perturbation Rate (s)",
           y = "Prediction Error",
           color = "Method") +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    
    plots[[i]] <- p
  }
  
  return(plots)
}

## Create plots
cat("\n=== Creating performance plots ===\n")
worst_plots <- create_performance_plot(worst_summary, "Worst Case")
ave_plots <- create_performance_plot(ave_summary, "Average Case")

## Save plots
cat("Saving plots...\n")

# Save worst case plots
pdf(file.path(output_dir, "worst_case_performance.pdf"), width = 10, height = 8)
for(p in worst_plots) {
  print(p)
}
dev.off()

# Save average case plots
pdf(file.path(output_dir, "average_case_performance.pdf"), width = 10, height = 8)
for(p in ave_plots) {
  print(p)
}
dev.off()

## Create comparison heatmaps
cat("\n=== Creating comparison heatmaps ===\n")

create_heatmap <- function(summary_df, method = "Ours", case_type = "Worst Case") {
  plot_data <- summary_df %>%
    select(smax, `true-s`, paste0(method, "_mean")) %>%
    rename(error = paste0(method, "_mean"))
  
  p <- ggplot(plot_data, aes(x = `true-s`, y = smax, fill = error)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = median(plot_data$error, na.rm = TRUE)) +
    labs(title = sprintf("%s - %s Method", case_type, method),
         x = "True Perturbation Rate (s)",
         y = "Fitted smax",
         fill = "Error") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# Create heatmaps for each method
methods <- c("Ours", "SS", "SA", "RA", "MI", "PA")

pdf(file.path(output_dir, "worst_case_heatmaps.pdf"), width = 8, height = 6)
for(method in methods) {
  p <- create_heatmap(worst_summary, method, "Worst Case")
  print(p)
}
dev.off()

pdf(file.path(output_dir, "average_case_heatmaps.pdf"), width = 8, height = 6)
for(method in methods) {
  p <- create_heatmap(ave_summary, method, "Average Case")
  print(p)
}
dev.off()

## Print key findings
cat("\n=== Key Findings ===\n")

# Find optimal smax for each true-s
optimal_smax_worst <- worst_summary %>%
  group_by(`true-s`) %>%
  slice_min(Ours_mean, n = 1) %>%
  select(`true-s`, smax, Ours_mean)

optimal_smax_ave <- ave_summary %>%
  group_by(`true-s`) %>%
  slice_min(Ours_mean, n = 1) %>%
  select(`true-s`, smax, Ours_mean)

cat("\nOptimal smax for each true-s (Worst Case):\n")
print(optimal_smax_worst)

cat("\nOptimal smax for each true-s (Average Case):\n")
print(optimal_smax_ave)

# Compare method performance at true-s = 0.5
cat("\n\nMethod comparison at high perturbation (true-s = 0.5):\n")
worst_high_pert <- worst_summary %>%
  filter(`true-s` == 0.5) %>%
  select(smax, Ours_mean, SS_mean, SA_mean, RA_mean, MI_mean, PA_mean)

cat("Worst Case:\n")
print(worst_high_pert)

ave_high_pert <- ave_summary %>%
  filter(`true-s` == 0.5) %>%
  select(smax, Ours_mean, SS_mean, SA_mean, RA_mean, MI_mean, PA_mean)

cat("\nAverage Case:\n")
print(ave_high_pert)

## Save session info
sink(file.path(output_dir, "session_info.txt"))
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Number of seeds analyzed:", worst_results$n_seeds, "\n\n")
cat("Parameters:\n")
print(worst_results$params)
cat("\n\nSession Info:\n")
print(sessionInfo())
sink()

cat("\n=== Analysis completed successfully ===\n")
cat("Results saved to:", output_dir, "\n")