############################################################
# File: Utils/Packages.R
# Purpose: Check and load required packages
# Usage: source("Utils/Packages.R")
############################################################

required_packages <- c(
  "data.table",
  "dplyr",
  "TwoSampleMR",
  "gwasvcf"
)

# Identify missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Load installed packages
invisible(lapply(setdiff(required_packages, missing_packages), function(p) {
  suppressPackageStartupMessages(
    library(p, character.only = TRUE)
  )
}))

# Report missing packages
if (length(missing_packages) > 0) {
  message("The following required packages are NOT installed:\n")
  message(paste0("  - ", missing_packages, collapse = "\n"))
  
  message("\nPlease install them before running the pipeline, e.g.:\n")
  message("install.packages(c(", 
          paste0('"', missing_packages, '"', collapse = ", "), 
          "))")
} else {
  message("All required packages are installed and loaded.")
}

# Return missing packages (invisibly)
invisible(missing_packages)