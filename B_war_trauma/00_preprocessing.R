################################################################################
#             PREPROCESSING - WAR/COMBAT TRAUMA (Axis B, GSE64813)
################################################################################
#namefile: 00_preprocessing.R
#creation dt: 18/08/2026
#last edit dt: 8/10/2026
#
# Axis B, step 1 of 3: import, sample ID mapping, time point filter,
# condition labels, low-expression filter, PCA inspection and normalization.
# Then run 01_DEG.R and 02_enrichment.R in the same R session.
#
# Dataset: GSE64813 (Breen et al., 2015), US Marines sampled before and after
# deployment; only the post-deployment samples are used (47 PTSD, 47 controls).
# The count matrix uses the sample descriptions (e.g. "Sample_1_1") as column
# names, so they are mapped to GEO IDs before the alignment.

# Imports & Libraries ----
source("00_config.R")
source("00_functions.R")
source("00_source.R")
load_pipeline_libraries()

cfg <- DIMENSIONS$war_trauma

plot_dir   <- file.path(cfg$output_dir, "plots")
data_dir   <- file.path(cfg$output_dir, "data")
tables_dir <- file.path(cfg$output_dir, "tables")
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Data import from GEO ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

raw <- load_geo_dataset(cfg$dataset_id, cfg$source_dir, cfg$raw_counts_pattern)

dim(raw$raw_counts) # 27974x188
head(colnames(raw$raw_counts))
head(rownames(raw$metadata))
length(intersect(colnames(raw$raw_counts), rownames(raw$metadata))) #188
head(raw$metadata[, c("title", "description")])

# Map sample descriptions to GEO IDs ----
desc_to_gsm <- setNames(rownames(raw$metadata), raw$metadata$description)
stopifnot(all(colnames(raw$raw_counts) %in% names(desc_to_gsm)))  # must be TRUE
colnames(raw$raw_counts) <- desc_to_gsm[colnames(raw$raw_counts)]
length(intersect(colnames(raw$raw_counts), rownames(raw$metadata))) # should be 188

# Time point filter: post-deployment samples only ----
if (!is.null(cfg$timepoint_filter)) {
  keep_gsm <- rownames(raw$metadata)[grepl(cfg$timepoint_filter, raw$metadata$title)]
  keep_samples <- colnames(raw$raw_counts) %in% keep_gsm
  table(keep_samples)
  raw$raw_counts <- raw$raw_counts[, keep_samples]
}
dim(raw$raw_counts) #27974x94

# Condition labels (Yes/No) ----
metadata <- build_condition_variable(raw$metadata,
                                      characteristics_col = cfg$characteristics_col,
                                      positive_values = cfg$positive_values,
                                      id_col = "geo_accession")
table(metadata$condition_cat)
table(metadata$condition_bin) # 94 vs 94: metadata not filtered by time point
head(rownames(metadata))
length(intersect(colnames(raw$raw_counts), rownames(metadata)))  # should be 94

# Sample alignment and low-expression filter ----
pre <- preprocess_axis_data(raw$raw_counts, metadata, plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)

dim(pre$counts)   # genes after filter:20956x94
dim(pre$meta)     # 94
table(pre$meta$condition_bin) # 47 vs. 47

# Library size by group (reported in Chapter 4): Wilcoxon rank-sum test
lib_M <- colSums(pre$counts) / 1e6                       # millions of counts
grp   <- pre$meta[colnames(pre$counts), cfg$condition_col]
summary(lib_M)                                           # min / median / max
tapply(lib_M, grp, median)                               # median by group
wilcox.test(lib_M ~ grp)                                 # p-value reported in the thesis
# Not significant: library sizes do not differ between groups.

# PCA inspection (no correction applied) ----
pca <- inspect_pca(pre$counts, pre$meta,
                    condition_col = cfg$condition_col,
                    group_label = cfg$group_label,
                    plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)

table(pca$meta$batch_candidate)

# Normalization (DESeq2 / VST), no batch correction ----
# The PCA before correction shows no clusters, so batch_col = NULL.
norm <- apply_normalization_pipeline(pre$counts, pre$meta,
                                      condition_col = cfg$condition_col,
                                      group_label = cfg$group_label,
                                      plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                                      batch_col = NULL)
dim(norm$vst) # PCA after normalization: plots/B00_PCA_post.png

# Objects kept in memory: raw, metadata, pre, pca, norm, cfg and the output
# folders. Next: 01_DEG.R, in the same R session.
