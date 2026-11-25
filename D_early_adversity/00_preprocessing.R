################################################################################
#           PREPROCESSING - EARLY-LIFE ADVERSITY (Axis D, GSE70603)
################################################################################
#namefile: 00_preprocessing.R
#creation dt: 15/08/2026
#
# Axis D, step 1 of 3: import, checks on control probes and time points,
# time point filter, condition labels, sex and age balance, genetic sex
# check, PCA inspection and normalization. Then run 01_DEG.R and
# 02_enrichment.R in the same R session.
#
# Dataset: GSE70603 (Schwaiger et al., 2016), Agilent microarray of CD14+
# monocytes. Differences from the RNA-seq axes:
#   1. processed log2 intensities are read from the series matrix
#      (load_geo_microarray_dataset()), not from the raw Agilent files;
#   2. normalization and batch correction with limma (quantile) and ComBat
#      instead of DESeq2 / ComBat_seq;
#   3. three samples per person (45 min before, 45 and 180 min after a
#      stressor): only the pre-stress samples (TBS45) are used.

# Imports & Libraries ----
source("00_config.R")
source("00_functions.R")
source("00_source.R")
load_pipeline_libraries()

cfg <- DIMENSIONS$early_adversity

plot_dir   <- file.path(cfg$output_dir, "plots")
data_dir   <- file.path(cfg$output_dir, "data")
tables_dir <- file.path(cfg$output_dir, "tables")
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Data import from GEO (processed microarray data) ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

raw <- load_geo_microarray_dataset(cfg$dataset_id, cfg$source_dir)

# age:ch1 comes as a string (e.g. "39y"): converted to numeric once here
raw$metadata$age <- as.numeric(gsub("y$", "", raw$metadata[["age:ch1"]]))

dim(raw$exprs)          # probes x samples (180 = 60 people x 3 time points)
colnames(raw$exprs)[1:10]
colnames(raw$metadata)
head(raw$metadata$title) # e.g. "A_TBS45_1": group_timepoint_subject
table(raw$metadata$characteristics_ch1)  # group labels
summary(raw$metadata$age)  # no NA expected

## Agilent control probes ----
# Control probes (positive/negative controls, spike-ins) must not be
# analyzed as genes. Expected: none in the series matrix.
grep("control|ERCC|DarkCorner", rownames(raw$exprs), ignore.case = TRUE, value = TRUE)
# 1. row names
head(rownames(raw$exprs), 10)
# 2. broader pattern for Agilent control probes
grep("^\\(\\+\\)|^\\(-\\)|Corner|^ETG|^RC[0-9]|NegativeControl|ERCC|E1A",
     rownames(raw$exprs), value = TRUE)
# 3. platform annotation: is there a control type column?
gpl <- GEOquery::getGEO("GPL17077", destdir = "D_early_adversity/01_source/GSE70603")
colnames(GEOquery::Table(gpl))
table(GEOquery::Table(gpl)$CONTROL_TYPE)   # if the column exists
ctrl <- GEOquery::Table(gpl)$ID[GEOquery::Table(gpl)$CONTROL_TYPE != "FALSE"]
sum(rownames(raw$exprs) %in% ctrl)   # expected: 0

## Group sizes and ages by time point ----
md <- raw$metadata
md$tp  <- sub("^[AB]_(TBS45|TAS45|TAS180)_.*$", "\\1", md$title)   # time point
md$grp <- ifelse(grepl("early adversity", md$source_name_ch1), "EA", "Control")

# 1. samples per group and time point (30 expected in each cell)
table(md$grp, md$tp)

# 2. age by group and time point: min, median, mean, max
aggregate(age ~ grp + tp, data = md,
          FUN = function(a) c(min = min(a), median = median(a),
                              mean = round(mean(a), 1), max = max(a)))

# 3. number of people under 45 at the pre-stress time point
with(md[md$tp == "TBS45", ], table(grp, age < 45))

# 4. do the ages of the same group match across time points?
same_ages <- function(g, t1, t2)
  identical(sort(md$age[md$grp == g & md$tp == t1]),
            sort(md$age[md$grp == g & md$tp == t2]))
same_ages("Control", "TBS45", "TAS180")   # TRUE expected
same_ages("Control", "TBS45", "TAS45")    # FALSE expected (groups swapped at +45 min)
same_ages("Control", "TAS45", "TBS45") ; same_ages("EA", "TBS45", "TAS180")
identical(sort(md$age[md$grp == "Control" & md$tp == "TAS45"]),
          sort(md$age[md$grp == "EA"      & md$tp == "TBS45"]))  # TRUE expected

# Time point filter: pre-stress samples only ----
keep_samples <- grepl(cfg$timepoint_filter, raw$metadata$title)
raw$exprs    <- raw$exprs[, keep_samples]
raw$metadata <- raw$metadata[keep_samples, ]
dim(raw$exprs)   # 60 samples expected (30 + 30)

# Condition labels (early adversity vs control -> Yes/No) ----
# The group is in source_name_ch1, not in characteristics_ch1 (see 00_config.R).
metadata <- build_condition_variable(raw$metadata,
                                      characteristics_col = cfg$characteristics_col,
                                      positive_values = cfg$positive_values,
                                      id_col = "geo_accession")
table(metadata$condition_cat)
table(metadata$condition_bin)

# Sex and age balance and genetic sex check ----
#   (1) sex x group: table and chi-square without Yates (Fisher if expected < 5)
#   (2) age by group: descriptive statistics and Welch t-test
#   (3) recorded sex vs expression of XIST and Y genes (RPS4Y1, DDX3Y,
#       KDM5D, UTY, EIF1AY) on the Agilent probes (GPL17077 annotation).
#       Check only: no sample is excluded here.
# Saves tables/D00_balance_sex_age.csv, tables/D00_sex_check.csv and
# plots/D00_sex_check.png
metadata$sex <- factor(metadata[["gender:ch1"]])
tb_sex <- table(metadata$sex, metadata[[cfg$condition_col]])
addmargins(tb_sex)
chi_sex <- test_2x2(tb_sex); unlist(chi_sex)   # chi-square without Yates, Fisher if expected < 5
tapply(metadata$age, metadata[[cfg$condition_col]], function(a) c(mean = mean(a), sd = sd(a), median = median(a), min = min(a), max = max(a)))
t_age <- t.test(metadata$age ~ metadata[[cfg$condition_col]]); t_age
write.csv(data.frame(test = c(paste0("sex_x_group_", chi_sex$test), "age_Welch_t"),
                     statistic = c(chi_sex$statistic, unname(t_age$statistic)),
                     df = c(chi_sex$df, unname(t_age$parameter)),
                     p_value = c(chi_sex$p, t_age$p.value)),
          file.path(tables_dir, paste0(cfg$axis_prefix, "00_balance_sex_age.csv")), row.names = FALSE)

gpl_D <- GEOquery::getGEO("GPL17077", AnnotGPL = TRUE, destdir = cfg$source_dir)
tab_D <- GEOquery::Table(gpl_D)
sym_col_D <- intersect(c("Gene symbol", "GENE_SYMBOL", "Symbol"), colnames(tab_D))[1]
sex_genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")
sex_probes <- tab_D[tab_D[[sym_col_D]] %in% sex_genes, c("ID", sym_col_D)]
sex_probes <- sex_probes[as.character(sex_probes$ID) %in% rownames(raw$exprs), ]
# probes found for each gene
table(sex_probes[[sym_col_D]])
samp_D <- intersect(colnames(raw$exprs), rownames(metadata))
x_sex  <- raw$exprs[as.character(sex_probes$ID), samp_D, drop = FALSE]   # log2 intensity
is_xist <- sex_probes[[sym_col_D]] == "XIST"
sex_check <- data.frame(sample   = samp_D,
                        title    = metadata[samp_D, "title"],
                        recorded = as.character(metadata[samp_D, "sex"]),
                        XIST     = colMeans(x_sex[is_xist, , drop = FALSE]),
                        Y_mean   = colMeans(x_sex[!is_xist, , drop = FALSE]))
sex_check$predicted <- ifelse(sex_check$Y_mean > sex_check$XIST, "male", "female")
# recorded vs predicted sex, and inconsistent samples
table(recorded = sex_check$recorded, predicted = sex_check$predicted)
sex_check[sex_check$recorded != sex_check$predicted, ]
write.csv(sex_check, file.path(tables_dir, paste0(cfg$axis_prefix, "00_sex_check.csv")), row.names = FALSE)
png(file.path(plot_dir, paste0(cfg$axis_prefix, "00_sex_check.png")), width = 900, height = 800, res = 130)
plot(sex_check$XIST, sex_check$Y_mean, pch = 19,
     col = ifelse(sex_check$recorded == "male", "steelblue", "#C0392B"),
     xlab = "XIST (log2 intensity)", ylab = "Y-chromosome genes (mean log2 intensity)",
     main = "Recorded sex vs sex-specific gene expression")
legend("topright", legend = c("recorded male", "recorded female"),
       col = c("steelblue", "#C0392B"), pch = 19, bty = "n")
dev.off()

# Sample alignment (microarray: no count filter) ----
pre <- preprocess_axis_data_microarray(raw$exprs, metadata, plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)
# Saves plots/D00_intensity_boxplot.png
dim(pre$exprs)
dim(pre$meta)

# Counterpart of the library size for arrays (reported in Chapter 4): median
# log2 intensity of each array before quantile normalization, compared
# between groups (range, median, Wilcoxon rank-sum test).
med_int <- apply(pre$exprs, 2, median)
grp     <- pre$meta[colnames(pre$exprs), cfg$condition_col]
summary(med_int)
tapply(med_int, grp, median)
wilcox.test(med_int ~ grp)

# PCA inspection (quantile normalization, no batch correction) ----
# batch_candidate (PC1 > 10) is only a heuristic suggestion, not used automatically.
pca <- inspect_pca_microarray(pre$exprs, pre$meta,
                               condition_col = cfg$condition_col,
                               group_label = cfg$group_label,
                               plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)

table(pca$meta$batch_candidate)

# Demographic variables on the PCA ----
# Sex (gender:ch1) and age (numeric, converted at import) tested against PC1.
pre$meta$sex <- factor(pre$meta[["gender:ch1"]])

pca_sex <- inspect_pca_microarray(pre$exprs, pre$meta,
                                   condition_col = cfg$condition_col,
                                   group_label = cfg$group_label,
                                   plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                                   extra_intgroup = "sex")

summary(aov(pca$pca$x[, 1] ~ pre$meta$sex))            # effect of sex on PC1
summary(lm(pca$pca$x[, 1] ~ pre$meta$age))              # effect of age on PC1

# Diagnostic PCA without X- and Y-chromosome probes ----
# Only to check whether other structure (technical batch, group) remains
# once sex is removed. The DEG analysis does not change: in limma, sex is
# a covariate (01_DEG.R).
hgnc_d <- read.delim("Transferability/Dictionary/hgnc_complete_set.txt", quote = "",
                     colClasses = "character")
xy_sym <- hgnc_d$symbol[grepl("^(X|Y)", hgnc_d$location)]
probe_sym <- tab_D[[sym_col_D]][match(rownames(pre$exprs), as.character(tab_D$ID))]
loc <- tab_D$CHROMOSOMAL_LOCATION[match(rownames(pre$exprs), as.character(tab_D$ID))]
# X/Y probes: HGNC symbol on X/Y or probe location on chrX/chrY (the
# latter also covers probes without a symbol; about 1,960 probes expected)
is_xy <- (!is.na(probe_sym) & probe_sym %in% xy_sym) | grepl("chr(X|Y):", loc)
sum(is_xy)                                   # X/Y probes left out of the plot
pca_noXY <- inspect_pca_microarray(pre$exprs[!is_xy, ], pre$meta,
                                    condition_col = cfg$condition_col,
                                    group_label = cfg$group_label,
                                    plot_dir = NULL, axis_prefix = cfg$axis_prefix)
summary(pca_noXY$pca)$importance[2, 1:2]                     # variance explained by PC1, PC2
summary(aov(pca_noXY$pca$x[, 1] ~ pre$meta$sex))             # sex no longer expected to dominate
summary(aov(pca_noXY$pca$x[, 1] ~ pre$meta[[cfg$condition_col]]))
summary(lm(pca_noXY$pca$x[, 1] ~ pre$meta$age))
save_pca_plot_generic(pca_noXY$pca, pre$meta[[cfg$condition_col]], cfg$group_label,
                      title = "PCA without X/Y-chromosome probes",
                      path = file.path(plot_dir, paste0(cfg$axis_prefix, "00_PCA_noXY.png")))

# Batch correction (ComBat) and normalization ----
# Quantile normalization is already done in inspect_pca_microarray().
# batch_col = NULL: no batch correction, decided after the PCA above.
norm <- apply_normalization_pipeline_microarray(pca$exprs_norm, pca$meta,
                                                 condition_col = cfg$condition_col,
                                                 group_label = cfg$group_label,
                                                 plot_dir = plot_dir,
                                                 axis_prefix = cfg$axis_prefix,
                                                 batch_col = NULL)
dim(norm$exprs)
unlist(chi_sex)


# Objects kept in memory: raw, metadata, pre, pca, norm, cfg and the output
# folders. Next: 01_DEG.R, in the same R session.
