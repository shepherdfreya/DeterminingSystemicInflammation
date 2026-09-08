############################################################
# File: Utils/Config.R
# Purpose: Project-wide constants, defaults, and helper settings
############################################################
# ---- PLINK ----
PLINK_BIN <- "/Users/lb21131/Documents/Programmes/plink/plink"

# ---- LD reference panel parent directory ----
AUX_DIR <- "Data/Aux"

# ---- Build-specific LD reference panel prefixes ----
# These should point to the PLINK prefix, not just the folder
REF_PANEL_PREFIXES <- list(
  GRCh37 = file.path(AUX_DIR, "EUR_HG19GRCh37", "EUR_HG19GRCh37"),
  GRCh38 = file.path(AUX_DIR, "EUR_HG38GRCh38", "EUR_HG38GRCh38")
)

# ---- Exposure preparation thresholds ----
EXPOSURE_MAF_THRESHOLD  <- 0.01
EXPOSURE_PVAL_THRESHOLD <- 5e-8
EXPOSURE_F_THRESHOLD    <- 10
EXPOSURE_CLUMP_R2       <- 0.001
EXPOSURE_CLUMP_KB       <- 10000

# ---- Proxy search settings ----
PROXY_TAG_R2   <- 0.8
PROXY_TAG_KB   <- 5000
PROXY_TAG_NSNP <- 5000

# ---- Colocalization settings ----
COLOC_H4_THRESHOLD <- 0.8
COLOC_MIN_OVERLAP  <- 10

# ---- PWCoCo colocalization settings ----
# PWCoCo binary and LD reference panels should be edited to match the local/HPC install.
PWCOCO_BIN <- Sys.getenv("PWCOCO_BIN", unset = "pwcoco")
PWCOCO_WINDOW <- 5e5
PWCOCO_MAF <- 0.01
PWCOCO_H4_THRESHOLD <- 0.8
PWCOCO_PPH4_SCALE <- 100
PWCOCO_REF_PANEL_PREFIXES <- REF_PANEL_PREFIXES

# Build GRCh37 MHC caution interval. Adjust if using a different convention.
MHC_CHR <- "6"
MHC_START <- 25000000
MHC_END <- 34000000

# ---- Evidence appraisal settings ----
EVIDENCE_ALPHA <- 0.05
EVIDENCE_BONFERRONI_N_TESTS <- NA



# ---- Cis vs Trans analysis ----
CIS_WINDOW <- 1e6
cis_window <- CIS_WINDOW
#Protein locations should be provided by user in trait_manifest

# Output subdirectories within the chosen output directory
SUBDIR_TABLES  <- "Tables"
SUBDIR_FIGURES <- "Figures"
SUBDIR_LOGS    <- "Logs"

# Plot settings
PALETTE_MAIN <- c(
  "#FF0000",
  "#00FF00",
  "#0000FF",
  "#FFFF00",
  "#00FFFF",
  "#FF00FF",
  "#FFA500",
  "#800080",
  "#A52A2A",
  "#808080"
)

theme_freya <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

# Global options
options(
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE
)
