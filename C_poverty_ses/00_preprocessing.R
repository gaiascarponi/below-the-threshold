################################################################################
#             PREPROCESSING - CHRONIC POVERTY (Axis C, GSE221766)
################################################################################
#namefile: 00_preprocessing.R
#creation dt: 15/08/2026
#
# Axis C, step 1 of 3: import, genetic sex check, time point filter,
# condition labels, low-expression filter, PCA inspection, covariate checks,
# removal of sex-chromosome genes and normalization. Then run 01_DEG.R and
# 02_enrichment.R in the same R session.
#
# Dataset: GSE221766 (HANDLS, Arnold et al.), RNA-seq of PBMCs. Two
# differences from Axis A:
#   1. one count file per sample, loaded with load_geo_dataset_persample();
#   2. two time points per person: only the first (T1) is used, so that the
#      samples are independent. T2 is used only for the sex check.

# Imports & Libraries ----
source("00_config.R")
source("00_functions.R")
source("00_source.R")
load_pipeline_libraries()

cfg <- DIMENSIONS$poverty_ses

plot_dir   <- file.path(cfg$output_dir, "plots")
data_dir   <- file.path(cfg$output_dir, "data")
tables_dir <- file.path(cfg$output_dir, "tables")
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Data import from GEO (one count file per sample) ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
raw <- load_geo_dataset_persample(cfg$dataset_id, cfg$source_dir, cfg$raw_counts_pattern)

dim(raw$raw_counts)     # genes x samples (408 = 204 people x 2 time points)
colnames(raw$raw_counts)[1:10]
colnames(raw$metadata)
table(raw$metadata[["poverty status:ch1"]])
sum(table(raw$metadata[["poverty status:ch1"]]))

# Genetic sex check on all 408 samples (T1 and T2) ----
# XIST (expressed only with two X chromosomes) and 5 Y-chromosome genes, in
# CPM. Done before the T1 filter: the T2 sample of the same person shows
# whether the error is in the T1 sample or in the recorded sex.
sex_ens <- c(XIST = "ENSG00000229807", RPS4Y1 = "ENSG00000129824",
             DDX3Y = "ENSG00000067048", KDM5D = "ENSG00000012817",
             UTY = "ENSG00000183878", EIF1AY = "ENSG00000198692")
cnt <- raw$raw_counts[!grepl("^__", rownames(raw$raw_counts)), ]   # without the HTSeq summary rows
ids <- sub("\\.[0-9]+$", "", rownames(cnt))
stopifnot(all(sex_ens %in% ids))                                  # all 6 genes must be present
cpm <- t(t(as.matrix(cnt[match(sex_ens, ids), ])) / colSums(cnt) * 1e6)
rownames(cpm) <- names(sex_ens)

# metadata rows are GSM IDs, matrix columns are sample titles (HA001-T1)
meta_idx  <- match(colnames(cnt), raw$metadata$title)
stopifnot(!any(is.na(meta_idx)))
sex_check <- data.frame(sample   = colnames(cnt),
                        gsm      = rownames(raw$metadata)[meta_idx],
                        XIST     = round(cpm["XIST", ], 1),
                        Y_total  = round(colSums(cpm[-1, , drop = FALSE]), 1),
                        recorded = raw$metadata[meta_idx, "Sex:ch1"])
# Thresholds: in consistent samples XIST is < 50 CPM in men and > 1000 in
# women, Y_total > 490 in men and < 55 in women, so 100 separates them well
sex_check$genetic <- with(sex_check,
  ifelse(XIST > 100 & Y_total > 100, "mixed",
  ifelse(XIST > 100, "Female",
  ifelse(Y_total > 100, "Male", "undetermined"))))
sex_check$consistent <- sex_check$genetic == sex_check$recorded
table(sex_check$recorded, sex_check$genetic)
sex_check[!sex_check$consistent, ]            # all inconsistent samples (T1 and T2)

# T1 samples inconsistent while their T2 is consistent: the error is in the T1 sample
bad_T1 <- sex_check$sample[!sex_check$consistent & grepl("-T1$", sex_check$sample)]
T2_ok  <- sex_check$consistent[match(sub("-T1$", "-T2", bad_T1), sex_check$sample)]
sex_mismatch <- bad_T1[!is.na(T2_ok) & T2_ok]
sex_mismatch                                  # expected: HA029-T1, HA139-T1, HA150-T1
write.csv(sex_check, file.path(tables_dir, paste0(cfg$axis_prefix, "00_sex_check.csv")), row.names = FALSE)
rm(cnt, ids, cpm, meta_idx)

# Time point filter: T1 only ----
if (!is.null(cfg$timepoint_filter)) {
  keep_samples <- grepl(cfg$timepoint_filter, colnames(raw$raw_counts))
  raw$raw_counts <- raw$raw_counts[, keep_samples]
}
dim(raw$raw_counts)   # 204 samples

table(raw$metadata[["cell type:ch1"]])  # constant (PBMC): not a confounder
table(raw$metadata[["tissue:ch1"]])     # constant (whole blood): not a confounder

# Condition labels (below vs above the poverty line -> Yes/No) ----
metadata <- build_condition_variable(raw$metadata,
                                      characteristics_col = cfg$characteristics_col,
                                      positive_values = cfg$positive_values,
                                      id_col = "title")
table(metadata$condition_cat)
table(metadata$condition_bin)

# Exclusion of 3 T1 samples with inconsistent genetic sex ----
# At T1, the time point analyzed here:
#   - HA139-T1 (GSM6895708): recorded male, purely female profile (high XIST,
#     no Y expression); its T2 is male, as recorded
#   - HA150-T1 (GSM6895730): recorded female, purely male profile; its T2 is
#     female, as recorded
#     -> in both cases the T1 sample does not belong to the recorded person
#        (sample swap or labeling error)
#   - HA029-T1 (GSM6895526): recorded male, mixed profile (XIST at about half
#     the female level and Y genes at about half the male level); its T2 is
#     purely male -> the T1 sample is probably contaminated with RNA from
#     another person (an XXY karyotype would also show at T2)
# All three are below the poverty line. sex_mismatch comes from the sex
# check above; the script stops if it differs from these 3 samples.
stopifnot(setequal(sex_mismatch, c("HA029-T1", "HA139-T1", "HA150-T1")))
sum(rownames(metadata) %in% sex_mismatch)
sum(colnames(raw$raw_counts) %in% sex_mismatch)
metadata <- metadata[!rownames(metadata) %in% sex_mismatch, ]
raw$raw_counts <- raw$raw_counts[, !colnames(raw$raw_counts) %in% sex_mismatch]
dim(raw$raw_counts)            # 201 samples expected
table(metadata$condition_bin[rownames(metadata) %in% colnames(raw$raw_counts)])   # 94 Yes / 107 No

# Quick check of the low-expression thresholds (the plots are saved by
# preprocess_axis_data() below)
dim(raw$raw_counts)

gene_sums <- rowSums(raw$raw_counts)
summary(gene_sums)
for (thr in c(1, 5, 10, 20, 50, 100)) {
  cat("Threshold >=", thr, "-> genes remaining:", sum(gene_sums >= thr),
      "out of", length(gene_sums), "\n")
}

lib_sizes <- colSums(raw$raw_counts)
summary(lib_sizes)

# Sample alignment and low-expression filter ----
# Also saves plots/C00_low_expression_filter.png and C00_library_size.png
pre <- preprocess_axis_data(raw$raw_counts, metadata, min_count_sum = 10, plot_dir = plot_dir,
                             axis_prefix = cfg$axis_prefix)

dim(pre$counts)
dim(pre$meta)

# Library size by group (reported in Chapter 4): range, median and
# Wilcoxon rank-sum test on the counts after the low-expression filter.
lib_M <- colSums(pre$counts) / 1e6                       # millions of counts
grp   <- pre$meta[colnames(pre$counts), cfg$condition_col]
summary(lib_M)                                           # min / median / max
tapply(lib_M, grp, median)                               # median by group
wilcox.test(lib_M ~ grp)                                 # p-value reported in the thesis


# PCA inspection (no correction applied) ----
# batch_candidate (PC1 > 10) is only a heuristic suggestion, not used automatically.
pca <- inspect_pca(pre$counts, pre$meta,
                    condition_col = cfg$condition_col,
                    group_label = cfg$group_label,
                    plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)

table(pca$meta$batch_candidate)

# Demographic and technical variables on the PCA ----
# Sex and race (categorical) and pct_assigned (continuous technical proxy
# from the HTSeq summary rows, see extract_htseq_qc_stats()).
pre$meta$sex  <- factor(pre$meta[["Sex:ch1"]])
pre$meta$race <- factor(pre$meta[["race:ch1"]])
pre$meta$pct_assigned <- raw$qc_stats[rownames(pre$meta), "pct_assigned"]
summary(pre$meta$pct_assigned)

pca_sex <- inspect_pca(pre$counts, pre$meta,
                        condition_col = cfg$condition_col,
                        group_label = cfg$group_label,
                        plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                        extra_intgroup = "sex")

pca_race <- inspect_pca(pre$counts, pre$meta,
                         condition_col = cfg$condition_col,
                         group_label = cfg$group_label,
                         plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                         extra_intgroup = "race")

# pct_assigned is continuous: plotted with a continuous color scale on the
# PCA already computed

save_pca_plot_continuous(pca$pcaData, pre$meta$pct_assigned, "% reads assigned",
                          title = "PCA by % reads assigned (technical QC proxy)",
                          path = file.path(plot_dir, paste0(cfg$axis_prefix, "00_PCA_pre_by_pct_assigned.png")))

summary(aov(pca$pcaData$PC1 ~ pre$meta$sex))            # p < 2e-16, R^2 about 96% (201 samples): sex dominates PC1
summary(aov(pca$pcaData$PC1 ~ pre$meta$race))            # p = 0.173 (201 samples): not significant
summary(lm(pca$pcaData$PC1 ~ pre$meta$pct_assigned))     # p = 0.123 (201 samples): not significant

# Sex and race balance between groups (to rule out confounding before using
# sex as a covariate)
table(pre$meta$sex, pre$meta[[cfg$condition_col]])
table(pre$meta$race, pre$meta[[cfg$condition_col]])
# independence tests below (test_2x2)

## Comparison runs (not part of the final pipeline, kept for reference) ----
# ComBat_seq batch correction with the PC1 heuristic, and with sex treated as
# a batch instead of a covariate. Results (FDR < 0.10):
#   - sex as covariate (data not changed)         -> 56-57 genes
#   - PC1 heuristic batch correction (ComBat_seq) -> 28 genes
#   - ComBat_seq on sex                           -> 22 genes
# Treating sex as a covariate keeps more power, so it is the choice used below.
# pre$meta$batch_pc1 <- pca$meta$batch_candidate
#
# norm_pc1 <- apply_normalization_pipeline(pre$counts, pre$meta,
#                                           condition_col = cfg$condition_col,
#                                           group_label = cfg$group_label,
#                                           plot_dir = plot_dir,
#                                           axis_prefix = paste0(cfg$axis_prefix, "_pc1test"),
#                                           batch_col = "batch_pc1")
#
# deg_pc1 <- run_deg_analysis(norm_pc1$dds, norm_pc1$vst, condition_col = cfg$condition_col,
#                              positive_level = "Yes", negative_level = "No",
#                              thr_pval = 0.10, thr_fc = 1)
#
# table(deg_pc1$res_filtered$diffexpressed)              # down=15, up=13 -> 28 genes
# sum(deg_pc1$res_filtered$padj < 0.10, na.rm = TRUE)
#
# norm_sexbatch <- apply_normalization_pipeline(pre$counts, pre$meta,
#                                                condition_col = cfg$condition_col,
#                                                group_label = cfg$group_label,
#                                                plot_dir = plot_dir,
#                                                axis_prefix = paste0(cfg$axis_prefix, "_sexbatchtest"),
#                                                batch_col = "sex")
#
# deg_sexbatch <- run_deg_analysis(norm_sexbatch$dds, norm_sexbatch$vst, condition_col = cfg$condition_col,
#                                   positive_level = "Yes", negative_level = "No",
#                                   thr_pval = 0.10, thr_fc = 1)
#
# table(deg_sexbatch$res_filtered$diffexpressed)         # down=10, up=12 -> 22 genes
# sum(deg_sexbatch$res_filtered$padj < 0.10, na.rm = TRUE)

# Sex and race balance by poverty status ----
# Tables and Pearson chi-square without Yates correction (Fisher's exact test
# if an expected count is below 5, see test_2x2()). If sex or race were
# unbalanced they would be confounders; if balanced, sex enters the design
# only as a covariate. Saves tables/C00_balance_chisq.csv
bal_vars <- list(sex_x_poverty  = c("sex",  cfg$condition_col),
                 race_x_poverty = c("race", cfg$condition_col),
                 sex_x_race     = c("sex",  "race"))
bal_tab <- do.call(rbind, lapply(names(bal_vars), function(nm) {
  v  <- bal_vars[[nm]]
  tb <- table(pre$meta[[v[1]]], pre$meta[[v[2]]])
  ct <- test_2x2(tb)
  print(nm); print(addmargins(tb)); print(unlist(ct))
  data.frame(comparison = nm,
             cells = paste(paste0(rep(rownames(tb), ncol(tb)), "/", rep(colnames(tb), each = nrow(tb)),
                                  "=", as.vector(tb)), collapse = "; "),
             test = ct$test, min_expected = ct$min_expected,
             chisq = ct$statistic, df = ct$df, p_value = ct$p)
}))
# Expected without Yates: sex and race vs poverty p about 0.92, sex x race
# p about 0.43; all chi-square, minimum expected counts > 45
bal_tab
write.csv(bal_tab, file.path(tables_dir, paste0(cfg$axis_prefix, "00_balance_chisq.csv")), row.names = FALSE)

# Removal of X- and Y-chromosome genes ----
# As in the original study (Arnold et al.), genes on the sex chromosomes are
# removed. With the design ~ sex + condition, the poverty log2FC of genes
# expressed in one sex only (XIST, Y genes) is an artifact: in an earlier run,
# 8 of 13 DEGs were XIST and Y genes. All X and Y genes are removed,
# pseudoautosomal ones included. Chromosome from HGNC (location) and NCBI
# gene_info (chromosome), via Ensembl ID. Done after the PCA, which documents
# the effect of sex, and before normalization.
hgnc_c <- read.delim("Transferability/Dictionary/hgnc_complete_set.txt", quote = "",
                     colClasses = "character")
ncbi_c <- read.delim(gzfile("Transferability/Dictionary/Homo_sapiens.gene_info.gz"),
                     quote = "", comment.char = "", colClasses = "character")
has_ens  <- grepl("ENSG[0-9]+", ncbi_c$dbXrefs)
ens_ncbi <- regmatches(ncbi_c$dbXrefs, regexpr("ENSG[0-9]+", ncbi_c$dbXrefs))
ncbi_xy  <- ens_ncbi[ncbi_c$chromosome[has_ens] %in% c("X", "Y", "X|Y")]
hgnc_xy  <- hgnc_c$ensembl_gene_id[grepl("^(X|Y)", hgnc_c$location) & hgnc_c$ensembl_gene_id != ""]
sexchr_ens <- unique(c(hgnc_xy, ncbi_xy))
ids_c <- sub("\\.[0-9]+$", "", rownames(pre$counts))
drop_sexspec <- ids_c %in% sexchr_ens
# Number of genes removed (a few hundred, about 3-4%), XIST and Y genes included
sum(drop_sexspec); round(100 * mean(drop_sexspec), 1)
c("ENSG00000229807", "ENSG00000012817", "ENSG00000129824", "ENSG00000067048") %in% ids_c[drop_sexspec]   # XIST, KDM5D, RPS4Y1, DDX3Y: TRUE expected
pre$counts <- pre$counts[!drop_sexspec, ]
dim(pre$counts)

# Normalization (DESeq2 / VST) with sex as a covariate ----
# No real technical batch was found (submission date and instrument are
# constant, pct_assigned is not associated with PC1), so no ComBat_seq
# (batch_col = NULL). Sex explains most of the PC1 variance (ANOVA,
# p < 2e-16) and is balanced between groups, so it enters the design as a
# covariate, not as a batch (see the comparison runs above).
norm <- apply_normalization_pipeline(pre$counts, pre$meta,
                                      condition_col = cfg$condition_col,
                                      group_label = cfg$group_label,
                                      plot_dir = plot_dir,
                                      axis_prefix = cfg$axis_prefix,
                                      batch_col = NULL,
                                      extra_covariates = "sex")

dim(norm$vst)
norm$dds@design

# Objects kept in memory: raw, metadata, pre, pca, norm, cfg and the output
# folders. Next: 01_DEG.R, in the same R session.
