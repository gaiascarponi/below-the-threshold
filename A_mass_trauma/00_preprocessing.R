################################################################################
#                  PREPROCESSING - MASS TRAUMA (Axis A, WTC)
################################################################################
#namefile: 00_preprocessing.R
#creation dt: 18/08/2026
#
# Axis A, step 1 of 3: import, condition labels, low-expression filter,
# batch correction and normalization. Then run 01_DEG.R and 02_enrichment.R
# in the same R session (they reuse pre, norm, cfg and the output folders).
# Axis A uses run_normalization_pipeline(), which estimates the batch from
# PC1 (PC1 > 10 on the VST scale) and corrects it with ComBat_seq. The
# other axes use apply_normalization_pipeline(), after inspecting the PCA
# with inspect_pca() on Axes C and D.

# Imports & Libraries ----
source("00_config.R")
source("00_functions.R")
source("00_source.R")
load_pipeline_libraries()
cfg <- DIMENSIONS$mass_trauma
plot_dir   <- file.path(cfg$output_dir, "plots")
data_dir   <- file.path(cfg$output_dir, "data")
tables_dir <- file.path(cfg$output_dir, "tables")
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Data import from GEO ----
# Downloads the data if they are not already in the source folder.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
raw <- load_geo_dataset(cfg$dataset_id, cfg$source_dir, cfg$raw_counts_pattern)
dim(raw$raw_counts)
head(raw$metadata)

# Condition labels (Yes/No) ----
metadata <- build_condition_variable(raw$metadata,
                                      characteristics_col = cfg$characteristics_col,
                                      positive_values = cfg$positive_values)
table(metadata$condition_cat)
table(metadata$condition_bin)

# Sample alignment and low-expression filter ----
pre <- preprocess_axis_data(raw$raw_counts, metadata, plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)
# Saves plots/A00_low_expression_filter.png and A00_library_size.png
dim(pre$counts)   # genes left after the filter
dim(pre$meta)     # aligned samples
# Library size by group (reported in Chapter 4): range, median and
lib_M <- colSums(pre$counts) / 1e6                       # millions of counts
grp   <- pre$meta[colnames(pre$counts), cfg$condition_col]
summary(lib_M)                                           # min / median / max
tapply(lib_M, grp, median)                               # median by group
wilcox.test(lib_M ~ grp)                                 # p-value reported in the thesis
# Wilcoxon rank-sum test on the counts after the low-expression filter.

# Batch correction (ComBat_seq) and normalization (DESeq2 / VST) ----
norm <- run_normalization_pipeline(pre$counts, pre$meta,
                                    condition_col = cfg$condition_col,
                                    group_label = cfg$group_label,
                                    plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)
# The batch is estimated from the PC1 heuristic and always corrected.
table(norm$meta$batch)
dim(norm$vst)

# Objects kept in memory: raw, metadata, pre, norm, cfg, plot_dir,
# data_dir, tables_dir. Next: 01_DEG.R, in the same R session.