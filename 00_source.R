################################################################################
#                         SOURCE CODE
################################################################################
#namefile: 00_source.R
#creation dt: 11/08/2026

#dir structure
folders <- c(
  "A_mass_trauma",
  "B_war_trauma",
  "C_poverty_ses",
  "D_early_adversity"
)
subfolders <- c(
  "01_source",      # raw data
  "02_processed",   # cleaned/normalized data
  "03_results",     # DEG tables, Hub genes, Plots...
  "04_network"      # Regulation and network analysis (SWIMmeR, MIENTURNET, STRING)
)

for (fold in folders) {
  for (sub in subfolders) {
    full_path <- file.path(fold, sub)
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = TRUE)
      message(paste("Folder created:", full_path))
    }
  }
}



