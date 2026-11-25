################################################################################
#                     BELOW THE THRESHOLD: FULL PIPELINE
################################################################################
# Run from the project root (open Project.Rproj), so that all paths are relative.
# Input data must be downloaded first (see README.md, "Data not included").
# Each script loads 00_config.R and 00_functions.R (or DALY/00_daly_functions.R)
# by itself, so every block can also be run on its own.

# MOLECULAR ANALYSIS ----
## LEVEL 1 ----

### Axis A ----
source("A_mass_trauma/00_preprocessing.R", echo = TRUE)
source("A_mass_trauma/01_DEG.R", echo = TRUE)
source("A_mass_trauma/02_enrichment.R", echo = TRUE)
### Axis B ----
source("B_war_trauma/00_preprocessing.R", echo = TRUE)
source("B_war_trauma/01_DEG.R", echo = TRUE)
source("B_war_trauma/02_enrichment.R", echo = TRUE)
### Axis C ----
source("C_poverty_ses/00_preprocessing.R", echo = TRUE)
source("C_poverty_ses/01_DEG.R", echo = TRUE)
source("C_poverty_ses/02_enrichment.R", echo = TRUE)
### Axis D ----
source("D_early_adversity/00_preprocessing.R", echo = TRUE)
source("D_early_adversity/01_DEG.R", echo = TRUE)
source("D_early_adversity/02_enrichment.R", echo = TRUE)

## LEVEL 2 ----
#Not run from here:
#  -SWIMmeR (external tool, seed 18, k = 3): see A_mass_trauma/04_network/swimmer/README.md
#  -MIENTURNET and STRING (web tools), on the gene lists in A_mass_trauma/04_network/input/

##LEVEL 3 ----

source("Transferability/00_harmonize_gene_ids.R", echo = TRUE)
source("Transferability/01_signature_transfer_gsea.R", echo = TRUE)
source("Transferability/02_functional_gsea_all_axes.R", echo = TRUE)
source("Transferability/plot_hallmark_heatmap.R", echo = TRUE)
source("Transferability/03_leading_edge.R", echo = TRUE)
source("Transferability/04_cell_composition.R", echo = TRUE)
source("Transferability/05_specificity.R", echo = TRUE)

# SOCIAL COST ANALYSIS ----
# DALY model; does not use the expression data. GBD 2023 inputs must be in
# DALY/data/source/ (see README.md). 00_daly_functions.R is loaded by each script.

source("DALY/01_axisA_PTSD_DALY.R", echo = TRUE)
source("DALY/02_axisB_PTSD_DALY.R", echo = TRUE)
source("DALY/03_axisC_poverty_DALY.R", echo = TRUE)
source("DALY/04_axisD_ACE_DALY.R", echo = TRUE)
source("DALY/05_compare_axes.R", echo = TRUE)
source("DALY/06_italy_comparison.R", echo = TRUE)
