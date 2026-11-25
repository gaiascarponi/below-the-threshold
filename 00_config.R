################################################################################
#                         CONFIG CODE
################################################################################
#namefile: 00_config.R
#creation dt: 11/08/2026

# Thresholds ----
# Default thresholds for all axes
THR_PVAL <- 0.05
THR_FC <- 1.1
PRC_IQR <- 0.2

# Datasets ----
# One entry per axis, with everything run_dimension_pipeline() needs: where the raw data are, how the binary condition is built and where
# the outputs go.
# Raw data are downloaded from GEO with GEOquery into
# <axis folder>/01_source/<dataset_id>/ if they are not already there:
#   - series matrix (metadata)       -> getGEO(dataset_id, ...)
#   - raw counts (supplementary file) -> getGEOSuppFiles(dataset_id, ...)
# raw_counts_pattern is the regex that picks the right supplementary file
# ("RawCounts" is usually enough).
#
# positive_values: values of characteristics_col that become "Yes" in
# condition_bin. Everything else becomes "No".

DIMENSIONS <- list(

  # Axis A: acute mass trauma (WTC responders)
  mass_trauma = list(
    name                = "mass_trauma",
    axis_prefix         = "A",  # prefix of the output file names (A_..., A00_...)
    dataset_id          = "GSE97356",
    source_dir          = "A_mass_trauma/01_source/GSE97356",
    raw_counts_pattern  = "RawCounts",
    characteristics_col = "characteristics_ch1",
    positive_values     = c("Current", "Past"),
    condition_col       = "condition_bin",
    group_label         = "PTSD",
    output_dir          = "A_mass_trauma/03_results",
    swimmer_dir         = "SWIMmeR/project/PTSD_WTC/dataset/PTSD",
    thr_pval            = THR_PVAL,
    thr_fc              = THR_FC,
    prc_iqr             = PRC_IQR
  ),

  # Axis C: chronic poverty (HANDLS)
  poverty_ses = list(
    name                = "poverty_ses",
    axis_prefix         = "C",
    dataset_id          = "GSE221766",
    source_dir          = "C_poverty_ses/01_source/GSE221766",
    raw_counts_pattern  = "\\.cnts\\.txt\\.gz$",  # one file per sample
    characteristics_col = "poverty status:ch1",   
    positive_values     = c("Below"),             # below the poverty line = Yes
    condition_col       = "condition_bin",
    group_label         = "Poverty",
    output_dir          = "C_poverty_ses/03_results",
    swimmer_dir         = "SWIMmeR/project/poverty_ses/dataset",
    timepoint_filter    = "-T1$",  # first time point only (the second is 4-5 years later)
    # Same thresholds as the other axes (the original study used FDR < 0.10).
    thr_pval            = THR_PVAL,
    thr_fc              = THR_FC,
    prc_iqr             = PRC_IQR
  ),

  # Axis B: war and combat trauma (US Marines)
  war_trauma = list(
    name                = "war_trauma",
    axis_prefix         = "B",
    dataset_id          = "GSE64813",
    source_dir          = "B_war_trauma/01_source/GSE64813",
    raw_counts_pattern  = "RAW_COUNT_MATRIX",
    characteristics_col = "condition:ch1",
    # The "Case" group has two labels depending on the time point:
    # "Case (PTSD risk)" before deployment, "Case (PTSD)" after it.
    positive_values     = c("Case (PTSD risk)", "Case (PTSD)"),
    condition_col       = "condition_bin",
    group_label         = "PTSD",
    output_dir          = "B_war_trauma/03_results",
    swimmer_dir         = "SWIMmeR/project/war_trauma/dataset",
    # One time point per axis, as for A and C: post-deployment, when the
    # PTSD status is known (94 subjects, 47 cases and 47 controls). Sample
    # titles look like "Case Post [Sample_1_3]", so the filter matches "Post".
    timepoint_filter    = "Post",
    thr_pval            = THR_PVAL,
    thr_fc              = THR_FC,
    prc_iqr             = PRC_IQR
  ),

  # Axis D: early-life adversity (CD14+ monocytes)
  early_adversity = list(
    name                = "early_adversity",
    axis_prefix         = "D",
    dataset_id          = "GSE70603",
    source_dir          = "D_early_adversity/01_source/GSE70603",
    platform            = "microarray",  # Agilent GPL17077: limma branch, no DESeq2 / ComBat_seq
    # The group is not in characteristics_ch1 (only tissue, age and sex) but
    # in source_name_ch1 (checked: 30 "control group", 30 "early adversity").
    characteristics_col = "source_name_ch1",
    positive_values     = c("Isolated monocytes, 45min pre-stress, early adversity"),
    condition_col       = "condition_bin",
    group_label         = "EarlyAdversity",
    output_dir          = "D_early_adversity/03_results",
    # Main analysis on the pre-stress time point (TBS45, 45 min before the
    # Trier Social Stress Test): a baseline group comparison, as on the other
    # axes. positive_values above only matches TBS45. The post-stress samples
    # (TAS45, TAS180) are not used here.
    timepoint_filter    = "_TBS45_",  # 45 min before the stressor
    swimmer_dir         = "SWIMmeR/project/early_adversity/dataset",
    thr_pval            = THR_PVAL,
    thr_fc              = THR_FC,
    prc_iqr             = PRC_IQR
  )
)