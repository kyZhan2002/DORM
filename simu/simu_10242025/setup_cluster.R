
## One-time setup script for cluster
## Run this before submitting simulation jobs to install all required packages

cat("=== Setting up renv environment on cluster ===\n")

# Load here package (should be available in base R or pre-installed)
if (!require("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

library(here)

# Find and activate renv
renv_activate_path <- here("renv", "activate.R")
if(file.exists(renv_activate_path)) {
  cat("Activating renv from:", renv_activate_path, "\n")
  source(renv_activate_path)
} else {
  stop("Error: renv/activate.R not found at ", renv_activate_path)
}

# Restore all packages from lockfile
cat("Restoring packages from renv.lock...\n")
renv::restore(prompt = FALSE)

cat("\n=== Setup completed successfully! ===\n")
cat("You can now submit simulation jobs.\n")