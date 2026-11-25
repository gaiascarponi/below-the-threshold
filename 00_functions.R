################################################################################
#                         FUNCTIONS CODE
################################################################################
#namefile: 00_functions.R
#creation dt: 11/08/2026

# Libraries ----
load_pipeline_libraries <- function() {
  suppressPackageStartupMessages({
    library(GEOquery)
    library(DESeq2)
    library(ggplot2)
    library(sva)
    library(pheatmap)
    library(ggpubr)
    library(dplyr)
    library(writexl)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(enrichR)
    library(stringr)
    library(limma) 
  })
}

# Download ----
# Download GEO supplementary files (RNA-seq or microarray) into dest_dir.
# Skips existing files and extracts any .tar archives.

download_geo_dataset_raw <- function(dataset_id, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  existing <- list.files(dest_dir)
  if (length(existing) == 0) {
    message("Downloading the supplementary files of ", dataset_id, " in ", dest_dir, " ...")
    getGEOSuppFiles(dataset_id, baseDir = dirname(dest_dir), makeDirectory = TRUE)
  } else {
    message(dataset_id, " already available locally (", length(existing), " file), skip download.")
  }

  # Extract any .tar archive (e.g. GSE221766_RAW.tar, GSE70603_RAW.tar)
  tar_files <- list.files(dest_dir, pattern = "\\.tar$", full.names = TRUE)
  for (tf in tar_files) {
    extract_dir <- file.path(dest_dir, sub("\\.tar$", "", basename(tf)))
    if (!dir.exists(extract_dir)) {
      dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
      message("Extracting ", basename(tf), " ...")
      untar(tf, exdir = extract_dir)
    }
  }

  # Also download the series matrix (metadata), whatever the platform
  series_matrix_path <- file.path(dest_dir, paste0(dataset_id, "_series_matrix.txt.gz"))
  if (!file.exists(series_matrix_path)) {
    gse <- tryCatch(getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = dest_dir),
                     error = function(e) { message("Series matrix not available: ", e$message); NULL })
  }

  message("Contents of ", dest_dir, ":")
  print(list.files(dest_dir, recursive = TRUE))
  invisible(list.files(dest_dir, recursive = TRUE, full.names = TRUE))
}

# Import  ----
# For datasets with one count file per sample  (GSE221766: each GSM has "<GSM>_<title>.cnts.txt.gz", two columns, gene
# and count, no header, inside the series RAW.tar). Merge the per-sample count files into one matrix (genes x samples).
# sample_name_fn takes the sample name from the file name (by default, everything before ".cnts.txt.gz")

merge_persample_counts <- function(dest_dir, file_pattern = "\\.cnts\\.txt\\.gz$",
                                    sample_name_fn = function(f) sub("\\.cnts\\.txt\\.gz$", "", basename(f))) {
  count_files <- list.files(dest_dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
  if (length(count_files) == 0) {
    stop("No per-sample count file found in ", dest_dir,
         " with pattern '", file_pattern, "'.")
  }

  counts_list <- lapply(count_files, function(f) {
    df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("gene", "count"))
    # HTSeq adds 5 alignment summary rows at the end of each file
    # (__no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique). They are not genes, so they are dropped.
    df <- df[!startsWith(df$gene, "__"), ]
    setNames(df$count, df$gene)
  })
  sample_names <- vapply(count_files, sample_name_fn, character(1))
  # GSM6895484_HA001-T1 -> keep the part after the first underscore (sample title)
  sample_names <- sub("^GSM[0-9]+_", "", sample_names)

  common_genes <- Reduce(intersect, lapply(counts_list, names))
  counts_mat <- sapply(counts_list, function(x) x[common_genes])
  rownames(counts_mat) <- common_genes
  colnames(counts_mat) <- sample_names
  return(as.data.frame(counts_mat))
}

# For each per-sample count file, read the 5 HTSeq summary rows dropped by merge_persample_counts() and compute a technical proxy per sample:
# pct_assigned, the share of reads assigned to a gene. It is not the metric of the original paper (% reads mapping to exons), but it is the only
# per-sample technical information available in the GEO files.
extract_htseq_qc_stats <- function(dest_dir, file_pattern = "\\.cnts\\.txt\\.gz$",
                                    sample_name_fn = function(f) sub("\\.cnts\\.txt\\.gz$", "", basename(f))) {
  count_files <- list.files(dest_dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
  if (length(count_files) == 0) {
    stop("No per-sample count file found in ", dest_dir,
         " with pattern '", file_pattern, "'.")
  }

  stats_list <- lapply(count_files, function(f) {
    df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("gene", "count"))
    is_special <- startsWith(df$gene, "__")
    special_sums <- setNames(df$count[is_special], df$gene[is_special])
    assigned <- sum(df$count[!is_special])
    total <- assigned + sum(special_sums)

    get_special <- function(name) if (name %in% names(special_sums)) special_sums[[name]] else NA_real_

    data.frame(
      assigned_reads       = assigned,
      no_feature            = get_special("__no_feature"),
      ambiguous              = get_special("__ambiguous"),
      too_low_aQual           = get_special("__too_low_aQual"),
      not_aligned              = get_special("__not_aligned"),
      alignment_not_unique      = get_special("__alignment_not_unique"),
      total_reads                = total,
      pct_assigned                = 100 * assigned / total
    )
  })

  sample_names <- vapply(count_files, sample_name_fn, character(1))
  sample_names <- sub("^GSM[0-9]+_", "", sample_names)

  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- sample_names
  return(stats_df)
}

# Like load_geo_dataset(), for datasets with per-sample counts: downloads the files if needed (download_geo_dataset_raw()) and merges them. 
# Also returns qc_stats, the per-sample proxy from extract_htseq_qc_stats(), indexed like colnames(raw_counts).
load_geo_dataset_persample <- function(dataset_id, source_dir, counts_pattern = "\\.cnts\\.txt\\.gz$") {
  download_geo_dataset_raw(dataset_id, source_dir)
  series_matrix_path <- file.path(source_dir, paste0(dataset_id, "_series_matrix.txt.gz"))
  gse <- getGEO(filename = series_matrix_path, getGPL = FALSE)
  metadata <- pData(gse)
  raw_counts <- merge_persample_counts(source_dir, counts_pattern)
  qc_stats <- extract_htseq_qc_stats(source_dir, counts_pattern)
  return(list(gse = gse, metadata = metadata, raw_counts = raw_counts, qc_stats = qc_stats))
}

# Import: microarray ----
# Axis D (GSE70603, Agilent). The processed log2 intensities are already in the series matrix (exprs(gse)), so the raw Agilent files are 
# not needed: only the series matrix is downloaded, as for the metadata of the other axes.

load_geo_microarray_dataset <- function(dataset_id, source_dir) {
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  series_matrix_path <- file.path(source_dir, paste0(dataset_id, "_series_matrix.txt.gz"))
  if (!file.exists(series_matrix_path)) {
    message("Series matrix not found locally, downloading it from GEO: ", dataset_id)
    gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = source_dir)
    gse <- gse[[1]]
  } else {
    gse <- getGEO(filename = series_matrix_path, getGPL = FALSE)
  }

  metadata <- pData(gse)
  exprs_mat <- exprs(gse)  # processed log2 intensities

  return(list(gse = gse, metadata = metadata, exprs = exprs_mat))
}

# Import and condition labels ----
# Load the series matrix (metadata) and the raw counts of a GEO dataset. Missing files are downloaded into source_dir:
#   - series matrix -> getGEO(dataset_id, destdir = source_dir)
#   - raw counts    -> getGEOSuppFiles(dataset_id, ...)
# raw_counts_pattern is a regex that picks the right supplementary file

load_geo_dataset <- function(dataset_id, source_dir, raw_counts_pattern = "RawCounts") {
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  # Metadata 
  series_matrix_path <- file.path(source_dir, paste0(dataset_id, "_series_matrix.txt.gz"))
  if (!file.exists(series_matrix_path)) {
    message("Series matrix not found locally, downloading it from GEO: ", dataset_id)
    gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = source_dir)
    gse <- gse[[1]]
  } else {
    gse <- getGEO(filename = series_matrix_path, getGPL = FALSE)
  }
  metadata <- pData(gse)
  # Raw Counts
  raw_counts_files <- list.files(source_dir, pattern = raw_counts_pattern, full.names = TRUE)
  if (length(raw_counts_files) == 0) {
    message("Raw counts not found locally, downloading them from GEO: ", dataset_id)
    getGEOSuppFiles(dataset_id, baseDir = dirname(source_dir), makeDirectory = TRUE,
                     filter_regex = raw_counts_pattern)
    raw_counts_files <- list.files(source_dir, pattern = raw_counts_pattern, full.names = TRUE)
  }
  if (length(raw_counts_files) == 0) {
    stop("No raw count file found or downloadable for ", dataset_id,
         " with pattern '", raw_counts_pattern, "'. Check the name of the supplementary file on GEO.")
  }
  # Detect the separator from the first line: some count files are comma-separated (Axes A and C), others tab-separated (Axis B,
  # GSE64813_RAW_COUNT_MATRIX.txt.gz). read.csv() on a tab-separated file reads each row as a single field and returns no data columns
  first_line <- readLines(raw_counts_files[1], n = 1)
  n_comma <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  n_tab   <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
  sep <- if (n_tab > n_comma) "\t" else ","
  message("load_geo_dataset(): separator detected for '", basename(raw_counts_files[1]),
          "': ", ifelse(sep == "\t", "tab", "comma"))
  raw_counts <- read.csv(raw_counts_files[1], row.names = 1, check.names = FALSE, sep = sep)
  return(list(gse = gse, metadata = metadata, raw_counts = raw_counts))
}

# Build the binary condition (e.g. PTSD Yes/No) from a text column of sample characteristics (characteristics_ch1 or similar)
build_condition_variable <- function(metadata,
                                      characteristics_col = "characteristics_ch1",
                                      positive_values,
                                      id_col = "title",
                                      strip_pattern = ".*: ") {
  metadata$condition_cat <- gsub(strip_pattern, "", metadata[[characteristics_col]])
  metadata$condition_bin <- factor(
    ifelse(metadata$condition_cat %in% positive_values, "Yes", "No"),
    levels = c("No", "Yes")
  )
  rownames(metadata) <- metadata[[id_col]]
  return(metadata)
}

# Preprocessing: sample alignment and low-expression filter ----
preprocess_axis_data <- function(raw_counts, metadata, min_count_sum = 10, plot_dir = NULL,
                                  axis_prefix = "") {
  common_samples <- intersect(colnames(raw_counts), rownames(metadata))
  raw_counts <- raw_counts[, common_samples]
  metadata <- metadata[common_samples, ]
  gene_sums <- rowSums(raw_counts)
  keep <- gene_sums >= min_count_sum
  # plot_dir is optional (default NULL: no plots). If given, two preprocessing diagnostics are saved for each RNA-seq axis:
  #   1. low_expression_filter.png: distribution of counts per gene with the filter threshold, on the aligned samples before the filter
  #   2. library_size.png: library size per sample, to spot samples with unusual sequencing depth
  if (!is.null(plot_dir)) {
    # Preprocessing plots get the "00" prefix (e.g. "A00_..."), so they sort apart from the other plots of the same axis
    png(file.path(plot_dir, paste0(axis_prefix, "00_low_expression_filter.png")), width = 800, height = 600, res = 100)
    hist(log10(gene_sums + 1), breaks = 100,
         main = "Distribution of total counts per gene (log10 scale)",
         xlab = "log10(summed raw count + 1)", col = "steelblue", border = "white")
    abline(v = log10(min_count_sum + 1), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste0("threshold = ", min_count_sum),
           col = "red", lty = 2, lwd = 2, bty = "n")
    dev.off()

    save_library_size_plot(raw_counts, file.path(plot_dir, paste0(axis_prefix, "00_library_size.png")))
  }

  return(list(counts = raw_counts[keep, ], meta = metadata))
}

# Preprocessing: microarray ----
# Sample alignment only. Log2 intensities are not counts, so there is no
# min_count_sum filter: uninformative probes are removed later by the same
# IQR filter used for RNA-seq (see run_deg_analysis_microarray()).

preprocess_axis_data_microarray <- function(exprs_mat, metadata, plot_dir = NULL, axis_prefix = "") {
  common_samples <- intersect(colnames(exprs_mat), rownames(metadata))
  exprs_mat <- exprs_mat[, common_samples]
  metadata <- metadata[common_samples, ]

  # Counterpart of library_size.png for arrays: the distribution of
  # intensities per array before quantile normalization, to spot outlier
  # arrays.
  if (!is.null(plot_dir)) {
    png(file.path(plot_dir, paste0(axis_prefix, "00_intensity_boxplot.png")), width = 900, height = 600, res = 100)
    boxplot(exprs_mat, main = "Distribution of raw log2 intensities per array",
            ylab = "log2 intensity", xlab = "Sample", las = 2, cex.axis = 0.5,
            col = "steelblue", outline = FALSE)
    dev.off()
  }

  return(list(exprs = exprs_mat, meta = metadata))
}

# Batch correction and normalization (VST) ----

# group_label is only used for the PCA plot labels (e.g. "PTSD")
run_normalization_pipeline <- function(counts, meta, condition_col = "condition_bin",
                                        group_label = condition_col,
                                        plot_dir = NULL, batch_pc1_thr = 10,
                                        axis_prefix = "") {

  dds_pre <- DESeqDataSetFromMatrix(countData = counts, colData = meta,
                                     design = as.formula(paste("~", condition_col)))
  vsd_pre <- vst(dds_pre, blind = FALSE)
  pcaData <- plotPCA(vsd_pre, intgroup = condition_col, returnData = TRUE)
  meta$batch <- ifelse(pcaData$PC1 > batch_pc1_thr, "Batch_A", "Batch_B")

  if (!is.null(plot_dir)) {
    save_pca_plot(vsd_pre, condition_col, group_label,
                   title = paste0("PCA of Gene Expression by ", group_label, " Status - Original Data"),
                   path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_pre.png")))
  }

  # ComBat_seq
  adj_counts <- ComBat_seq(as.matrix(counts), batch = meta$batch, group = meta[[condition_col]])

  # Final DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = adj_counts, colData = meta,
                                 design = as.formula(paste("~", condition_col)))
  dds <- DESeq(dds)
  vsd_post <- vst(dds, blind = FALSE)
  vst_data <- assay(vsd_post)

  if (!is.null(plot_dir)) {
    save_pca_plot(vsd_post, condition_col, group_label,
                   title = paste0("PCA of Gene Expression by ", group_label, " Status - Adjusted for Batch Effects"),
                   path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_post.png")))
    save_pca_plot(vsd_post, c(condition_col, "batch"), group_label,
                   title = "PCA: Condition and Batch (Post-Correction)",
                   path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_batch_check.png")))
  }

  return(list(dds = dds, vst = vst_data, meta = meta, vsd_post = vsd_post))
}

# PCA inspection ----
# Kept apart from the correction, so that the PCA can be inspected before
# deciding whether and how to correct for batch. The data are not changed:
# it computes a VST, draws the PCA colored by condition_col and, if asked,
# returns a heuristic batch assignment (PC1 > batch_pc1_thr) as an extra
# column of meta. That column is only a candidate: it is used only if it is
# passed to apply_normalization_pipeline().
inspect_pca <- function(counts, meta, condition_col = "condition_bin",
                         group_label = condition_col, plot_dir = NULL,
                         batch_pc1_thr = 10, extra_intgroup = NULL,
                         axis_prefix = "") {

  dds_pre <- DESeqDataSetFromMatrix(countData = counts, colData = meta,
                                     design = as.formula(paste("~", condition_col)))
  vsd_pre <- vst(dds_pre, blind = FALSE)

  pcaData <- plotPCA(vsd_pre, intgroup = condition_col, returnData = TRUE)
  meta$batch_candidate <- ifelse(pcaData$PC1 > batch_pc1_thr, "Batch_A", "Batch_B")

  if (!is.null(plot_dir)) {
    save_pca_plot(vsd_pre, condition_col, group_label,
                   title = paste0("PCA of Gene Expression by ", group_label, " Status - Before Correction"),
                   path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_pre.png")))

    # To color the PCA by another variable only (e.g. sex or race, not
    # combined with condition_col), pass its column name in extra_intgroup.
    if (!is.null(extra_intgroup)) {
      save_pca_plot(vsd_pre, extra_intgroup, extra_intgroup,
                     title = paste0("PCA by ", extra_intgroup),
                     path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_pre_by_", extra_intgroup, ".png")))
    }
  }

  return(list(vsd_pre = vsd_pre, pcaData = pcaData, meta = meta))
}

# Batch correction and normalization: explicit ----
# To be called after inspect_pca(). Unlike run_normalization_pipeline(),
# the batch is not estimated here: batch_col = NULL skips ComBat_seq
# (normalization only); otherwise it names the column of meta to use as
# batch (e.g. "batch_candidate" from inspect_pca()).
apply_normalization_pipeline <- function(counts, meta, condition_col = "condition_bin",
                                          group_label = condition_col,
                                          plot_dir = NULL, batch_col = NULL,
                                          axis_prefix = "", extra_covariates = NULL) {

  if (!is.null(batch_col)) {
    adj_counts <- ComBat_seq(as.matrix(counts), batch = meta[[batch_col]],
                              group = meta[[condition_col]])
  } else {
    adj_counts <- as.matrix(counts)
  }

  # extra_covariates: columns of meta added to the design before
  # condition_col (e.g. c("sex")), to adjust for confounders without
  # changing the counts (unlike ComBat_seq). condition_col stays last, so
  # the contrast tested in run_deg_analysis() does not change.
  design_terms <- c(extra_covariates, condition_col)
  design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

  dds <- DESeqDataSetFromMatrix(countData = adj_counts, colData = meta,
                                 design = design_formula)
  dds <- DESeq(dds)
  vsd_post <- vst(dds, blind = FALSE)
  vst_data <- assay(vsd_post)

  if (!is.null(plot_dir)) {
    correction_label <- if (is.null(batch_col)) "No Batch Correction" else "Adjusted for Batch Effects"
    save_pca_plot(vsd_post, condition_col, group_label,
                   title = paste0("PCA of Gene Expression by ", group_label, " Status - ", correction_label),
                   path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_post.png")))
    if (!is.null(batch_col)) {
      save_pca_plot(vsd_post, c(condition_col, batch_col), group_label,
                     title = "PCA: Condition and Batch (Post-Correction)",
                     path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_batch_check.png")))
    }
  }

  return(list(dds = dds, vst = vst_data, meta = meta, vsd_post = vsd_post))
}

# Batch correction and normalization: microarray ----
# Same logic as the RNA-seq branch (PCA to estimate the batch, correction,
# PCA check), with tools for continuous log2 data instead of counts:
#   - normalizeBetweenArrays() (limma) instead of VST
#   - ComBat() (sva) instead of ComBat_seq()

# top_var_rows: the most variable rows (probes) of a matrix.
# Mirrors the default of DESeq2's plotPCA() (ntop = 500), so that the PCA of
# Axis D is computed like the RNA-seq ones: the 500 most variable rows,
# centered and not scaled.
top_var_rows <- function(mat, ntop = 500) {
  rv <- apply(mat, 1, var)
  mat[order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))], , drop = FALSE]
}

run_normalization_pipeline_microarray <- function(exprs_mat, meta, condition_col = "condition_bin",
                                                    group_label = condition_col,
                                                    plot_dir = NULL, batch_pc1_thr = 10,
                                                    axis_prefix = "") {

  exprs_norm <- normalizeBetweenArrays(exprs_mat, method = "quantile")

  pca_pre <- prcomp(t(top_var_rows(exprs_norm)), scale. = FALSE)  # as in plotPCA(): 500 most variable probes, not scaled
  meta$batch <- ifelse(pca_pre$x[, 1] > batch_pc1_thr, "Batch_A", "Batch_B")

  if (!is.null(plot_dir)) {
    save_pca_plot_generic(pca_pre, meta[[condition_col]], group_label,
                           title = paste0("PCA of Gene Expression by ", group_label, " Status - Original Data"),
                           path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_pre.png")))
  }

  mod <- model.matrix(as.formula(paste("~", condition_col)), data = meta)
  exprs_adj <- ComBat(dat = exprs_norm, batch = meta$batch, mod = mod, par.prior = TRUE)

  if (!is.null(plot_dir)) {
    pca_post <- prcomp(t(top_var_rows(exprs_adj)), scale. = FALSE)  # as in plotPCA(): 500 most variable probes, not scaled
    save_pca_plot_generic(pca_post, meta[[condition_col]], group_label,
                           title = paste0("PCA of Gene Expression by ", group_label, " Status - Adjusted for Batch Effects"),
                           path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_post.png")))
  }

  return(list(exprs = exprs_adj, meta = meta))
}

# PCA inspection and explicit correction: microarray ----
# Same logic as inspect_pca() / apply_normalization_pipeline(), for
# microarray data: quantile normalization instead of VST, ComBat() instead
# of ComBat_seq(). run_normalization_pipeline_microarray() above corrects
# automatically from the PC1 heuristic, without inspection.

inspect_pca_microarray <- function(exprs_mat, meta, condition_col = "condition_bin",
                                    group_label = condition_col, plot_dir = NULL,
                                    batch_pc1_thr = 10, extra_intgroup = NULL,
                                    axis_prefix = "") {

  exprs_norm <- normalizeBetweenArrays(exprs_mat, method = "quantile")
  pca_pre <- prcomp(t(top_var_rows(exprs_norm)), scale. = FALSE)  # as in plotPCA(): 500 most variable probes, not scaled
  meta$batch_candidate <- ifelse(pca_pre$x[, 1] > batch_pc1_thr, "Batch_A", "Batch_B")

  if (!is.null(plot_dir)) {
    save_pca_plot_generic(pca_pre, meta[[condition_col]], group_label,
                           title = paste0("PCA of Gene Expression by ", group_label, " Status - Before Batch Correction"),
                           path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_pre.png")))

    # As in inspect_pca(): to color the PCA by another variable only (e.g.
    # sex or age), pass its column name in extra_intgroup.
    if (!is.null(extra_intgroup)) {
      save_pca_plot_generic(pca_pre, meta[[extra_intgroup]], extra_intgroup,
                             title = paste0("PCA by ", extra_intgroup),
                             path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_pre_by_", extra_intgroup, ".png")))
    }
  }

  return(list(exprs_norm = exprs_norm, pca = pca_pre, meta = meta))
}

# To be called after inspect_pca_microarray(). exprs_norm is the
# quantile-normalized matrix returned by inspect_pca_microarray(), not the
# raw one: quantile normalization is always applied, while the batch
# correction (ComBat) is optional.
apply_normalization_pipeline_microarray <- function(exprs_norm, meta, condition_col = "condition_bin",
                                                      group_label = condition_col,
                                                      plot_dir = NULL, batch_col = NULL,
                                                      axis_prefix = "") {

  if (!is.null(batch_col)) {
    mod <- model.matrix(as.formula(paste("~", condition_col)), data = meta)
    exprs_adj <- ComBat(dat = exprs_norm, batch = meta[[batch_col]], mod = mod, par.prior = TRUE)
  } else {
    exprs_adj <- exprs_norm
  }

  if (!is.null(plot_dir)) {
    pca_post <- prcomp(t(top_var_rows(exprs_adj)), scale. = FALSE)  # as in plotPCA(): 500 most variable probes, not scaled
    correction_label <- if (is.null(batch_col)) "No Batch Correction" else "Adjusted for Batch Effects"
    save_pca_plot_generic(pca_post, meta[[condition_col]], group_label,
                           title = paste0("PCA of Gene Expression by ", group_label, " Status - ", correction_label),
                           path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_post.png")))
    if (!is.null(batch_col)) {
      # Check after correction: PCA colored by the batch used.
      save_pca_plot_generic(pca_post, meta[[batch_col]], batch_col,
                             title = "PCA: Batch (Post-Correction)",
                             path = file.path(plot_dir, paste0(axis_prefix, "00_PCA_batch_check.png")))
    }
  }

  return(list(exprs = exprs_adj, meta = meta))
}

# Generic version of save_pca_plot() for any expression matrix (prcomp),
# used by the microarray branch, which has no DESeqTransform object.
save_pca_plot_generic <- function(pca, group_var, group_label, title, path,
                                   width = 800, height = 600, res = 100) {
  var_expl <- summary(pca)$importance[2, 1:2] * 100
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = group_var)
  png(path, width = width, height = height, res = res)
  p <- ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2) +
    labs(title = title,
         x = paste0("PC1: ", round(var_expl[1], 1), "% varianza"),
         y = paste0("PC2: ", round(var_expl[2], 1), "% varianza"),
         color = group_label)
  print(p)
  dev.off()
}

# DEG identification: IQR filter and DESeq2 results ----

run_deg_analysis <- function(dds, vst_data, condition_col = "condition_bin",
                              positive_level = "Yes", negative_level = "No",
                              prc_iqr = PRC_IQR, thr_fc = THR_FC, thr_pval = THR_PVAL,
                              plot_dir = NULL, axis_prefix = "", iqr_plot_xlim = NULL) {

  iqr_values <- apply(vst_data, 1, IQR)
  thr_perc <- quantile(iqr_values, prc_iqr)
  genes_pass_iqr <- names(iqr_values[iqr_values > thr_perc])

  # plot_dir is optional: if given, the histogram of the IQR distribution
  # is saved. iqr_plot_xlim (e.g. c(0, 6)) limits the x axis of that plot.
  if (!is.null(plot_dir)) {
    plot_iqr_distribution(iqr_values, thr_perc, prc_iqr, file.path(plot_dir, paste0(axis_prefix, "01_IQR_distribution.png")),
                           xlim = iqr_plot_xlim)
  }

  # alpha = thr_pval: DESeq2 uses alpha to tune its independent filtering
  # (genes with too low mean expression get padj = NA). The default is 0.1;
  # the DESeq2 documentation recommends matching it to the significance
  # threshold actually used.
  res_df <- as.data.frame(results(dds, contrast = c(condition_col, positive_level, negative_level),
                                   alpha = thr_pval))
  res_filtered <- res_df[rownames(res_df) %in% genes_pass_iqr, ]

  res_filtered$diffexpressed <- "delete"
  res_filtered$diffexpressed[res_filtered$padj <= thr_pval & res_filtered$log2FoldChange > log2(thr_fc)] <- "up"
  res_filtered$diffexpressed[res_filtered$padj <= thr_pval & res_filtered$log2FoldChange < -log2(thr_fc)] <- "down"

  sig_genes <- rownames(res_filtered[res_filtered$diffexpressed != "delete", ])

  export_tab <- res_filtered[res_filtered$diffexpressed != "delete", c("log2FoldChange", "padj")]
  export_tab$Gene_Name <- rownames(export_tab)
  export_tab <- export_tab[, c("Gene_Name", "padj", "log2FoldChange")]
  export_tab$cat <- ifelse(export_tab$log2FoldChange > 0, "UP", "DOWN")

  return(list(
    iqr_values = iqr_values, thr_perc = thr_perc,
    res_filtered = res_filtered, sig_genes = sig_genes,
    export_tab = export_tab
  ))
}

# Level 3: signature transfer (GSEA) ----
# Instead of requiring significant DEGs on every axis, test whether a
# signature defined on a source axis (Axis A) is enriched at the top or
# bottom of the complete gene ranking of a target axis (fgsea).

# Export the complete results table (all tested genes, unlike
# export_deg_table() / export_top10_xlsx()). It is the input of
# build_ranked_stats() when the transfer test runs in a separate R session.
export_deg_full_table <- function(res_filtered, path_rds) {
  saveRDS(res_filtered, path_rds)
}

# Map a results table indexed by Agilent probe ID (rownames) to gene
# symbols, with the platform annotation downloaded from GEO. Needed for
# Axis D (GPL17077) before build_ranked_stats() / run_go_enrichment().
# When several probes map to the same symbol, the one with the largest
# |stat| is kept (see below).
map_probes_to_symbols_agilent <- function(res_filtered, platform_id = "GPL17077") {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("The 'GEOquery' package is not installed.")
  }
  gpl <- GEOquery::getGEO(platform_id, AnnotGPL = TRUE, destdir = tempdir())
  tbl <- GEOquery::Table(gpl)

  symbol_col <- intersect(c("Gene symbol", "GENE_SYMBOL", "Symbol", "gene_symbol"), colnames(tbl))
  if (length(symbol_col) == 0) {
    stop("Symbol column not found in the annotation table of platform ", platform_id,
         ". Colonne disponibili: ", paste(colnames(tbl), collapse = ", "))
  }
  symbol_col <- symbol_col[1]
  id_col <- intersect(c("ID", "SPOT_ID"), colnames(tbl))[1]

  probe_to_symbol <- setNames(as.character(tbl[[symbol_col]]), as.character(tbl[[id_col]]))

  probe_ids <- rownames(res_filtered)
  symbols <- probe_to_symbol[probe_ids]
  # some platforms join several symbols with "///": keep the first
  symbols <- sub("///.*$", "", symbols)
  symbols <- trimws(symbols)

  keep <- !is.na(symbols) & symbols != ""
  res_mapped <- res_filtered[keep, , drop = FALSE]
  symbols_kept <- symbols[keep]

  # Several probes can map to the same symbol (common on Agilent arrays):
  # keep the probe with the largest |stat|, so that rownames are unique.
  stat_col <- intersect(c("stat", "t"), colnames(res_mapped))[1]
  ord <- order(-abs(res_mapped[[stat_col]]))
  res_mapped <- res_mapped[ord, , drop = FALSE]
  symbols_kept <- symbols_kept[ord]
  dup <- duplicated(symbols_kept)
  n_before <- nrow(res_mapped)
  res_mapped <- res_mapped[!dup, , drop = FALSE]
  symbols_kept <- symbols_kept[!dup]
  rownames(res_mapped) <- symbols_kept

  message(sprintf("Probes mapped to a gene symbol: %d / %d (%.1f%%); after deduplication (1 probe per symbol, max |stat|): %d unique symbols",
                   n_before, length(probe_ids), 100 * n_before / length(probe_ids), nrow(res_mapped)))
  return(res_mapped)
}

# Build the ranked vector for fgsea() from a results table (RNA-seq or
# microarray, both have a "stat" column): names are gene symbols, values
# are the test statistic (Wald for DESeq2, moderated t for limma), sorted in
# decreasing order. Ensembl IDs (Axis C, "ENSG" prefix) are mapped to
# symbols; symbols (Axis A) and mapped probes (Axis D) are kept as they are.
# Duplicates after mapping keep the entry with the largest |stat|.
build_ranked_stats <- function(res_filtered, stat_col = "stat") {
  # Harmonized table (harmonize_deg_table(), column "entrez"): ranking
  # indexed by Entrez Gene ID, one row per gene.
  if ("entrez" %in% colnames(res_filtered)) {
    ok <- !is.na(res_filtered$entrez) & !is.na(res_filtered[[stat_col]])
    ranked <- setNames(res_filtered[[stat_col]][ok], res_filtered$entrez[ok])
    ranked <- ranked[!duplicated(names(ranked))]
    return(sort(ranked, decreasing = TRUE))
  }
  ids <- rownames(res_filtered)
  is_ensembl <- mean(grepl("^ENSG", ids)) > 0.5

  if (is_ensembl) {
    gene_map <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL",
                             column = "SYMBOL", multiVals = "first")
    )
    symbols <- unname(gene_map)
  } else {
    symbols <- ids
  }

  stat_vals <- res_filtered[[stat_col]]
  df <- data.frame(symbol = symbols, stat = stat_vals, stringsAsFactors = FALSE)
  df <- df[!is.na(df$symbol) & !is.na(df$stat), ]
  df <- df[order(-abs(df$stat)), ]
  df <- df[!duplicated(df$symbol), ]

  ranked <- setNames(df$stat, df$symbol)
  ranked <- sort(ranked, decreasing = TRUE)
  return(ranked)
}

# Run the GSEA test (fgsea) for one or more signatures (e.g. "up",
# "down" and "all" from the Axis A DEGs) against one ranked vector (the
# target axis). Returns the standard fgsea table, one row per signature.
run_signature_transfer_gsea <- function(signatures, ranked_stats, min_size = 5, seed = 42) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("The 'fgsea' package is not installed. Install it with: BiocManager::install('fgsea')")
  }
  set.seed(seed)
  # nproc = 1: with parallel workers fgsea can give slightly different
  # results between runs even with set.seed(); with one process the result
  # is exactly reproducible.
  res <- fgsea::fgsea(pathways = signatures, stats = ranked_stats,
                       minSize = min_size, eps = 0, nproc = 1)
  res <- res[order(res$pval), ]
  return(as.data.frame(res))
}

# Gene identifier harmonization (Entrez Gene ID) ----
# The datasets use gene nomenclatures of different years (A: RefGene,
# around 2017; B: hg19/RSEM, around 2015; C: Ensembl; D: Agilent GPL17077
# annotation), so the same gene can appear under different symbols (e.g.
# CTSL vs CTSL1, KMT2D vs MLL2). For the cross-axis comparison, genes are
# matched through NCBI Entrez Gene IDs, which are stable over time.
# Dictionaries in Transferability/Dictionary/ (see its README):
#   - Homo_sapiens.gene_info.gz (NCBI Gene)
#   - hgnc_complete_set.txt     (HGNC)
# The conversion is applied after the DEG analysis, on the complete results
# tables, so the DEG results of each axis do not depend on gene names.

# Build the conversion tables from NCBI gene_info and the HGNC complete set.
build_gene_dictionary <- function(dict_dir = "Dictionary") {
  ncbi <- read.delim(gzfile(file.path(dict_dir, "Homo_sapiens.gene_info.gz")),
                     quote = "", comment.char = "", check.names = FALSE,
                     colClasses = "character")
  colnames(ncbi)[1] <- "tax_id"   # the original header is "#tax_id"
  hgnc <- read.delim(file.path(dict_dir, "hgnc_complete_set.txt"),
                     quote = "", check.names = FALSE, colClasses = "character")
  strip_q <- function(x) gsub('"', "", x, fixed = TRUE)

  # 1. official symbol -> Entrez (NCBI)
  sym2id <- setNames(ncbi$GeneID, ncbi$Symbol)
  sym2id <- sym2id[!duplicated(names(sym2id))]

  # 2. aliases and previous symbols -> Entrez (NCBI synonyms, HGNC alias/prev)
  split_pairs <- function(keys, ids) {
    lst <- strsplit(keys, "|", fixed = TRUE)
    data.frame(alias = unlist(lst), id = rep(ids, lengths(lst)), stringsAsFactors = FALSE)
  }
  al <- rbind(
    split_pairs(ncbi$Synonyms, ncbi$GeneID),
    split_pairs(strip_q(hgnc$alias_symbol), hgnc$entrez_id),
    split_pairs(strip_q(hgnc$prev_symbol),  hgnc$entrez_id)
  )
  al <- unique(al[al$alias != "" & al$alias != "-" & al$id != "", ])
  n_ids <- tapply(al$id, al$alias, function(x) length(unique(x)))
  alias_ambiguous <- names(n_ids)[n_ids > 1]          # alias -> several genes: dropped
  al_u <- al[al$alias %in% names(n_ids)[n_ids == 1], ]
  alias2id <- setNames(al_u$id, al_u$alias)
  alias2id <- alias2id[!duplicated(names(alias2id))]

  # 2b. HGNC previous symbols only -> Entrez (unique). Used to recover gene
  #     names corrupted by Excel (see to_entrez()): there the original name
  #     is the official symbol of its time (e.g. SEPT2), now an HGNC previous
  #     symbol (SEPTIN2). Some are ambiguous in the general alias table (SEPT2
  #     is also an NCBI alias of SEPTIN6), but not here.
  pv <- split_pairs(strip_q(hgnc$prev_symbol), hgnc$entrez_id)
  pv <- unique(pv[pv$alias != "" & pv$id != "", ])
  n_pv <- tapply(pv$id, pv$alias, function(x) length(unique(x)))
  pv <- pv[pv$alias %in% names(n_pv)[n_pv == 1], ]
  prev2id <- setNames(pv$id, pv$alias)

  # 3. Ensembl -> Entrez (NCBI dbXrefs, HGNC ensembl_gene_id), unique only
  xr <- strsplit(ncbi$dbXrefs, "|", fixed = TRUE)
  ens_ncbi <- data.frame(ens = unlist(xr), id = rep(ncbi$GeneID, lengths(xr)),
                         stringsAsFactors = FALSE)
  ens_ncbi <- ens_ncbi[startsWith(ens_ncbi$ens, "Ensembl:"), ]
  ens_ncbi$ens <- sub("^Ensembl:", "", ens_ncbi$ens)
  ens <- unique(rbind(ens_ncbi,
                      data.frame(ens = hgnc$ensembl_gene_id, id = hgnc$entrez_id,
                                 stringsAsFactors = FALSE)))
  ens <- ens[ens$ens != "" & ens$id != "", ]
  n_e <- tapply(ens$id, ens$ens, function(x) length(unique(x)))
  ens <- ens[ens$ens %in% names(n_e)[n_e == 1], ]
  ens2id <- setNames(ens$id, ens$ens)
  ens2id <- ens2id[!duplicated(names(ens2id))]

  # Entrez -> current symbol (for readable plots and leading edges)
  id2sym <- setNames(ncbi$Symbol, ncbi$GeneID)

  list(sym2id = sym2id, alias2id = alias2id, alias_ambiguous = alias_ambiguous,
       prev2id = prev2id,
       ens2id = ens2id, id2sym = id2sym,
       versions = c(ncbi_file = "Homo_sapiens.gene_info.gz", hgnc_file = "hgnc_complete_set.txt"))
}

# Convert a vector of identifiers (symbols or Ensembl IDs) to Entrez IDs.
# Returns a data.frame: input, entrez (NA if not convertible) and method
# (ensembl / symbol / alias / loc_id / excel_recovered / excel_corrupted / ambiguous / unmapped).
to_entrez <- function(ids, dict) {
  ids <- as.character(ids)
  out <- data.frame(input = ids, entrez = NA_character_, method = "unmapped",
                    stringsAsFactors = FALSE)
  excel_pat <- "^[0-9]+-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"

  is_ens <- grepl("^ENSG[0-9]+", ids)
  ens_clean <- sub("\\.[0-9]+$", "", ids)            # drop the version suffix .N
  hit <- is_ens & ens_clean %in% names(dict$ens2id)
  out$entrez[hit] <- dict$ens2id[ens_clean[hit]]; out$method[hit] <- "ensembl"

  todo <- !is_ens
  exc <- todo & grepl(excel_pat, ids, ignore.case = TRUE)
  out$method[exc] <- "excel_corrupted"; todo <- todo & !exc

  # Recover gene names corrupted by Excel (e.g. "2-Sep" from SEPT2): rebuild
  # the official symbol of the time and convert it with the HGNC previous
  # symbols (prev2id) or, failing that, with the unique aliases. Rules:
  #   N-Sep -> SEPT<N> (now SEPTIN<N>), except 15-Sep -> SEP15 (now SELENOF)
  #   N-Mar -> MARCH<N> (now MARCHF<N>) for N >= 3
  #   N-Dec -> DEC<N>   (1-Dec -> DEC1, now DELEC1)
  # 1-Mar and 2-Mar are not recovered: they can be MARCH1/2 (MARCHF1/2) or
  # MARC1/2 (MTARC1/2), so they stay "excel_corrupted" and are excluded.
  num <- suppressWarnings(as.integer(sub("-.*$", "", ids)))
  mon <- tolower(sub("^[0-9]+-", "", ids))
  legacy <- rep(NA_character_, length(ids))
  legacy[exc & mon == "sep"] <- ifelse(num[exc & mon == "sep"] == 15, "SEP15",
                                       paste0("SEPT", num[exc & mon == "sep"]))
  sel <- exc & mon == "mar" & !is.na(num) & num >= 3
  legacy[sel] <- paste0("MARCH", num[sel])
  sel <- exc & mon == "dec"
  legacy[sel] <- paste0("DEC", num[sel])
  rec_id <- ifelse(legacy %in% names(dict$prev2id), dict$prev2id[legacy],
            ifelse(legacy %in% names(dict$alias2id), dict$alias2id[legacy], NA_character_))
  hit <- exc & !is.na(rec_id)
  out$entrez[hit] <- rec_id[hit]; out$method[hit] <- "excel_recovered"

  hit <- todo & ids %in% names(dict$sym2id)
  out$entrez[hit] <- dict$sym2id[ids[hit]]; out$method[hit] <- "symbol"; todo <- todo & !hit

  hit <- todo & ids %in% names(dict$alias2id)
  out$entrez[hit] <- dict$alias2id[ids[hit]]; out$method[hit] <- "alias"; todo <- todo & !hit

  # "LOC<number>": the number is the Entrez Gene ID
  loc_num <- sub("^LOC([0-9]+)$", "\\1", ids)
  hit <- todo & grepl("^LOC[0-9]+$", ids) & loc_num %in% dict$sym2id
  out$entrez[hit] <- loc_num[hit]; out$method[hit] <- "loc_id"; todo <- todo & !hit

  amb <- todo & ids %in% dict$alias_ambiguous
  out$method[amb] <- "ambiguous"
  out
}

# Add entrez / symbol_current / map_method columns to a DEG results table
# (rownames = symbols or Ensembl IDs), drop the rows that cannot be
# converted and, when several rows map to the same Entrez ID, keep the one
# with the largest |stat| (same rule as build_ranked_stats()).
harmonize_deg_table <- function(res, dict, stat_col = "stat") {
  m <- to_entrez(rownames(res), dict)
  res$original_id    <- rownames(res)
  res$entrez         <- m$entrez
  res$map_method     <- m$method
  res$symbol_current <- unname(dict$id2sym[m$entrez])
  summary_tab <- table(factor(m$method, levels = c("ensembl", "symbol", "alias", "loc_id",
                                                   "excel_recovered", "excel_corrupted",
                                                   "ambiguous", "unmapped")))
  kept <- res[!is.na(res$entrez) & !is.na(res[[stat_col]]), ]
  kept <- kept[order(-abs(kept[[stat_col]])), ]
  n_dup <- sum(duplicated(kept$entrez))
  kept <- kept[!duplicated(kept$entrez), ]
  rownames(kept) <- kept$entrez
  attr(kept, "mapping_summary") <- summary_tab
  attr(kept, "n_collapsed_duplicates") <- n_dup
  kept
}

# DEG identification: microarray (limma) ----
# Same IQR filter and thresholds (THR_PVAL, THR_FC, PRC_IQR) as the RNA-seq
# branch, and the same output fields as run_deg_analysis() (res_filtered
# with log2FoldChange / padj / diffexpressed, iqr_values, thr_perc,
# sig_genes, export_tab), so the plotting, export and enrichment functions
# work unchanged.

run_deg_analysis_microarray <- function(exprs_adj, meta, condition_col = "condition_bin",
                                         positive_level = "Yes", negative_level = "No",
                                         prc_iqr = PRC_IQR, thr_fc = THR_FC, thr_pval = THR_PVAL,
                                         plot_dir = NULL, axis_prefix = "", iqr_plot_xlim = NULL,
                                         extra_covariates = NULL) {

  iqr_values <- apply(exprs_adj, 1, IQR)
  thr_perc <- quantile(iqr_values, prc_iqr)
  genes_pass_iqr <- names(iqr_values[iqr_values > thr_perc])

  if (!is.null(plot_dir)) {
    plot_iqr_distribution(iqr_values, thr_perc, prc_iqr, file.path(plot_dir, paste0(axis_prefix, "01_IQR_distribution.png")),
                           xlim = iqr_plot_xlim)   # as in run_deg_analysis()
  }

  # extra_covariates: as in apply_normalization_pipeline() for RNA-seq,
  # e.g. c("sex") -> design ~ sex + condition_col. condition_col stays last,
  # so the tested coefficient does not change.
  design_terms <- c(extra_covariates, condition_col)
  design <- model.matrix(as.formula(paste("~", paste(design_terms, collapse = " + "))), data = meta)
  fit <- lmFit(exprs_adj, design)
  fit <- eBayes(fit)

  # Coefficient of the positive level (condition_bin with levels
  # c("No", "Yes") -> column "condition_binYes" of the design matrix)
  coef_name <- paste0(condition_col, positive_level)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")

  # stat = limma moderated t (tt$t), the counterpart of the Wald statistic of
  # the RNA-seq branch (results(dds)$stat): both branches have a "stat"
  # column to rank genes for GSEA (see build_ranked_stats()).
  res_df <- data.frame(log2FoldChange = tt$logFC, pvalue = tt$P.Value, padj = tt$adj.P.Val, stat = tt$t,
                       row.names = rownames(tt))   # pvalue is used by the volcano plot and the histogram
  res_filtered <- res_df[rownames(res_df) %in% genes_pass_iqr, ]

  res_filtered$diffexpressed <- "delete"
  res_filtered$diffexpressed[res_filtered$padj <= thr_pval & res_filtered$log2FoldChange > log2(thr_fc)] <- "up"
  res_filtered$diffexpressed[res_filtered$padj <= thr_pval & res_filtered$log2FoldChange < -log2(thr_fc)] <- "down"

  sig_genes <- rownames(res_filtered[res_filtered$diffexpressed != "delete", ])

  export_tab <- res_filtered[res_filtered$diffexpressed != "delete", c("log2FoldChange", "padj")]
  export_tab$Gene_Name <- rownames(export_tab)
  export_tab <- export_tab[, c("Gene_Name", "padj", "log2FoldChange")]
  export_tab$cat <- ifelse(export_tab$log2FoldChange > 0, "UP", "DOWN")

  return(list(
    iqr_values = iqr_values, thr_perc = thr_perc,
    res_filtered = res_filtered, sig_genes = sig_genes,
    export_tab = export_tab
  ))
}

# Plots ----

save_pca_plot <- function(vsd, intgroup, group_label, title, path,
                           width = 800, height = 600, res = 100) {
  png(path, width = width, height = height, res = res)
  p <- plotPCA(vsd, intgroup = intgroup) + ggtitle(title)
  if (length(intgroup) == 1) p <- p + labs(color = group_label)
  print(p)
  dev.off()
}

# Like save_pca_plot(), for a continuous variable (e.g. % assigned reads).
# plotPCA() treats intgroup as a factor, so here pcaData (from
# plotPCA(..., returnData = TRUE), e.g. inspect_pca()$pcaData) is plotted
# directly with a continuous color scale.
save_pca_plot_continuous <- function(pcaData, values, var_label, title, path,
                                      width = 800, height = 600, res = 100) {
  df <- pcaData
  df$value <- values
  var_expl <- attr(pcaData, "percentVar")

  png(path, width = width, height = height, res = res)
  p <- ggplot(df, aes(x = PC1, y = PC2, color = value)) +
    geom_point(size = 2) +
    scale_color_gradient(low = "blue", high = "red") +
    ggtitle(title) +
    labs(color = var_label,
         x = if (!is.null(var_expl)) paste0("PC1: ", round(100 * var_expl[1]), "% variance") else "PC1",
         y = if (!is.null(var_expl)) paste0("PC2: ", round(100 * var_expl[2]), "% variance") else "PC2")
  print(p)
  dev.off()
}

# Library size (sum of raw counts) per sample, in increasing order, to spot
# samples with unusually low or high sequencing depth. Computed on the raw
# counts, before any gene is removed.
save_library_size_plot <- function(raw_counts, path, width = 900, height = 700, res = 120) {
  # Y axis in millions; sample names hidden (unreadable with hundreds of bars).
  lib_sizes <- colSums(raw_counts)
  png(path, width = width, height = height, res = res)
  barplot(sort(lib_sizes) / 1e6, main = "Library size per sample",
          ylab = "Total counts (millions)",
          xlab = "Samples (sorted by library size)",
          names.arg = NA, border = NA, col = "grey40")
  dev.off()
}

plot_iqr_distribution <- function(iqr_values, thr_perc, prc_iqr, path,
                                   width = 700, height = 800, res = 100, xlim = NULL) {
  png(path, width = width, height = height, res = res)
  # hist() sets xlim from the breaks only when xlim is not passed at all:
  # xlim = NULL is taken literally and fails in plot.window(), so the
  # arguments are built conditionally.
  hist_args <- list(x = iqr_values, breaks = 100,
                     main = "IQR distribution by gene",
                     xlab = "IQR value",
                     col = "steelblue",
                     border = "white")
  if (!is.null(xlim)) hist_args$xlim <- xlim
  do.call(hist, hist_args)
  abline(v = thr_perc, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste0(prc_iqr * 100, "th percentile"), col = "red", lty = 2, lwd = 2)
  dev.off()
}

save_volcano_plot <- function(res_df, path, thr_pval = THR_PVAL, thr_fc = THR_FC,
                               title = "Volcano plot", width = 800, height = 600, res = 100) {
  # The x axis follows the data (a fixed [-1.5, 1.5] would silently drop
  # genes with |log2FC| > 1.5), the y axis always shows the -log10(thr_pval)
  # line, even when no gene passes it, and the legend labels are readable.
  max_logFC <- max(abs(res_df$log2FoldChange), na.rm = TRUE)
  y_vals <- -log10(res_df$padj)
  max_y <- max(c(-log10(thr_pval), y_vals[is.finite(y_vals)]), na.rm = TRUE)
  png(path, width = width, height = height, res = res)
  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
    geom_point() +
    scale_x_continuous(limits = c(-max_logFC, max_logFC) * 1.05, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max_y * 1.1), expand = c(0, 0)) +
    scale_color_manual(values = c("delete" = "darkgrey", "up" = "yellow3", "down" = "blue"),   # yellow, as in the heatmap
                       labels = c("delete" = "Not significant", "up" = "Up", "down" = "Down")) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
          legend.title = element_text(colour = "black", size = 10, face = "bold"),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.box.background = element_rect(colour = "black")) +
    labs(title = title, x = "Log2 Fold-change (FC)", y = "-Log10(adj. p-value)", color = NULL) +
    geom_hline(yintercept = -log10(thr_pval), linetype = "dashed", color = "black") +
    geom_vline(xintercept = log2(thr_fc), linetype = "dashed", color = "black") +
    geom_vline(xintercept = -log2(thr_fc), linetype = "dashed", color = "black")
  print(volcano)
  dev.off()
}

save_heatmap_plot <- function(matrix_expr, path, annotation_col, anno_colors,
                               title = "Heatmap DEGs", show_rownames = NULL,
                               width = 800, height = 600, res = 100,
                               cluster_cols = TRUE, gaps_col = NULL) {
  # With 0 or 1 DEG the row-scaled heatmap cannot be drawn (a single gene
  # gives a vector without dim): skip it with a message.
  if (is.null(dim(matrix_expr)) || nrow(matrix_expr) < 2) {
    message("Heatmap not drawn: fewer than 2 significant genes (", path, ")")
    return(invisible(NULL))
  }
  if (is.null(show_rownames)) show_rownames <- nrow(matrix_expr) < 50
  png(path, width = width, height = height, res = res)
  pheatmap(matrix_expr,
           scale = "row",
           border_color = NA,
           cluster_cols = cluster_cols,
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "complete",
           annotation_col = annotation_col,
           annotation_colors = anno_colors,
           color = colorRampPalette(colors = c("blue", "blue3", "black", "yellow3", "yellow"))(100),
           show_rownames = show_rownames,
           show_colnames = FALSE,
           cutree_cols = if (cluster_cols) 2 else NA,
           gaps_col = gaps_col,
           cutree_rows = 2,
           main = title)
  dev.off()
}

save_boxplot_top_gene <- function(vst_data, res_filtered, meta, condition_col, group_label,
                                   direction = c("up", "down"), path,
                                   width = 800, height = 600, res = 100) {
  direction <- match.arg(direction)
  sub_res <- if (direction == "up") {
    res_filtered[res_filtered$log2FoldChange > 0, ]
  } else {
    res_filtered[res_filtered$log2FoldChange < 0, ]
  }
  top_gene <- rownames(sub_res)[which.min(sub_res$padj)]
  df_box <- data.frame(Expression = vst_data[top_gene, ], Group = meta[[condition_col]])
  title_txt <- if (direction == "up") paste("Most significant gene:", top_gene) else
    paste("Most significant DOWN-regulated gene:", top_gene)

  png(path, width = width, height = height, res = res)
  p <- ggplot(df_box, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 1) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = title_txt, x = group_label)
  print(p)
  dev.off()
  return(top_gene)
}

# Export: DEG tables and SWIMmeR input ----

export_deg_table <- function(export_tab, path_txt) {
  write.table(export_tab[, c("Gene_Name", "padj", "log2FoldChange")], file = path_txt,
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Gene lists (one symbol per line, ordered by padj) for MIENTURNET and
# STRING (Level 2). Writes <prefix>_DEG_all.txt, <prefix>_DEG_up.txt and
# <prefix>_DEG_down.txt in out_dir. With 0 DEGs (Axes B, D) no file is
# written and a message is shown.
export_deg_lists <- function(export_tab, out_dir, prefix) {
  if (is.null(export_tab) || nrow(export_tab) == 0) {
    message("export_deg_lists(", prefix, "): no DEGs, lists not written.")
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- export_tab[order(export_tab$padj), ]
  up   <- tab$Gene_Name[tab$log2FoldChange > 0]
  down <- tab$Gene_Name[tab$log2FoldChange < 0]
  writeLines(tab$Gene_Name, file.path(out_dir, paste0(prefix, "_DEG_all.txt")))
  writeLines(up,   file.path(out_dir, paste0(prefix, "_DEG_up.txt")))
  writeLines(down, file.path(out_dir, paste0(prefix, "_DEG_down.txt")))
  message("export_deg_lists(", prefix, "): ", length(up) + length(down), " DEG (",
          length(up), " up, ", length(down), " down) -> ", out_dir)
  invisible(list(all = tab$Gene_Name, up = up, down = down))
}

export_top10_xlsx <- function(export_tab, path_xlsx) {
  # dplyr:: is explicit because AnnotationDbi::select() (loaded with
  # clusterProfiler / org.Hs.eg.db) masks dplyr::select() and does not
  # handle a plain data.frame; same risk for filter() and arrange().
  top_10_up <- export_tab %>%
    dplyr::filter(cat == "UP") %>%
    dplyr::arrange(padj) %>%
    head(10) %>%
    dplyr::select(Gene_Name, log2FoldChange, padj)
  top_10_down <- export_tab %>%
    dplyr::filter(cat == "DOWN") %>%
    dplyr::arrange(padj) %>%
    head(10) %>%
    dplyr::select(Gene_Name, log2FoldChange, padj)
  write_xlsx(list("Top_10_UP" = top_10_up, "Top_10_DOWN" = top_10_down), path = path_xlsx)
}

# SWIMmeR (SWItchMiner, Paci et al., 2017) finds switch genes in a
# co-expression network (Pearson correlation) comparing two conditions. It
# does not need repeated measures and has been used on cross-sectional
# case-control studies, so it can be applied to any axis.
#
# expr_matrix is a normalized, batch-corrected expression matrix on the
# LINEAR scale (genes x samples):
#   - RNA-seq (A, B, C): counts(dds, normalized = TRUE), already linear
#   - microarray (D): 2^(norm$exprs), since limma/ComBat data are log2
export_swimmer_files <- function(meta, condition_col, expr_matrix, swimmer_dir, positive_level = "Yes", negative_level = "No") {
  dir.create(file.path(swimmer_dir, "list"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(swimmer_dir, "matrix"), recursive = TRUE, showWarnings = FALSE)

  samples_case <- rownames(meta[meta[[condition_col]] == positive_level, ])
  samples_ctrl <- rownames(meta[meta[[condition_col]] == negative_level, ])
  write.table(samples_case, file = file.path(swimmer_dir, "list", paste0(condition_col, "_", positive_level, "_list.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(samples_ctrl, file = file.path(swimmer_dir, "list", paste0(condition_col, "_", negative_level, "_list.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  swim_matrix <- as.data.frame(expr_matrix)
  swim_matrix <- cbind(Gene = rownames(swim_matrix), swim_matrix)
  write.table(swim_matrix, file = file.path(swimmer_dir, "matrix", "Normalized_Linear_Matrix_SWIM.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Functional enrichment (GO / enrichR) ----

run_go_enrichment <- function(genes, direction = c("up", "down"), plot_dir,
                               gene_list_fc = NULL, label = "", axis_prefix = "",
                               universe = NULL) {
  direction <- match.arg(direction)
  tag <- toupper(direction)
  dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "WikiPathways_2019_Human")

  # With no input genes (e.g. Axis B, no significant gene)
  # mean(grepl(...)) is NaN and the is_ensembl check fails: skip the
  # enrichment, as done below when length(symbols) == 0.
  if (length(genes) == 0) {
    warning(paste0("run_go_enrichment(", tag, "): no input genes, ",
                    "enrichment skipped for this direction."))
    return(list(enrichR = NULL, go_bp = NULL, go_cc = NULL, go_mf = NULL))
  }

  # Genes can come as Ensembl IDs (Axis C, "ENSG...") or as symbols
  # (Axis A, e.g. "FKBP5"). The format is detected from the share of
  # rownames starting with "ENSG": Ensembl IDs are mapped to symbols (needed
  # by Enrichr and enrichGO(keyType = "SYMBOL")), symbols are used as they
  # are (mapIds(keytype = "ENSEMBL") on symbols would return only NA).
  is_ensembl <- mean(grepl("^ENSG", genes)) > 0.5

  if (is_ensembl) {
    gene_map <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL",
                             column = "SYMBOL", multiVals = "first")
    )
    symbols <- unique(unname(gene_map[!is.na(gene_map)]))
    n_unmapped <- sum(is.na(gene_map))
    if (n_unmapped > 0) {
      warning(paste0("run_go_enrichment(", tag, "): ", n_unmapped, " genes out of ",
                      length(genes), " without a matching SYMBOL (excluded)."))
    }
  } else {
    symbols <- unique(genes)
  }

  if (length(symbols) == 0) {
    warning(paste0("run_go_enrichment(", tag, "): no gene mapped to SYMBOL, ",
                    "enrichment skipped for this direction."))
    return(list(enrichR = NULL, go_bp = NULL, go_cc = NULL, go_mf = NULL))
  }

  # universe (optional): the genes actually tested (e.g. rownames of
  # deg$res_filtered with non-NA padj), used as the enrichGO background
  # instead of all annotated genes. Without it, genes highly expressed in
  # blood (mitochondrial, ribosomal, immune) look enriched just because they
  # are the ones that can be measured. Default NULL keeps the whole-genome
  # background. Enrichr (web) does not use this background.
  universe_symbols <- NULL
  if (!is.null(universe) && length(universe) > 0) {
    if (mean(grepl("^ENSG", universe)) > 0.5) {
      u_map <- suppressMessages(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = universe, keytype = "ENSEMBL",
                               column = "SYMBOL", multiVals = "first")
      )
      universe_symbols <- unique(unname(u_map[!is.na(u_map)]))
    } else {
      universe_symbols <- unique(universe)
    }
    # input genes must be part of the background
    universe_symbols <- union(universe_symbols, symbols)
    message("run_go_enrichment(", tag, "): background (universe) = ",
            length(universe_symbols), " tested genes.")
  }
  go_args <- list(gene = symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  pAdjustMethod = "BH")
  if (!is.null(universe_symbols)) go_args$universe <- universe_symbols

  # gene_list_fc (for the cnetplot) uses the same IDs as `genes`: map it to
  # symbols only when needed (Ensembl).
  if (!is.null(gene_list_fc)) {
    if (is_ensembl) {
      fc_map <- suppressMessages(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(gene_list_fc), keytype = "ENSEMBL",
                               column = "SYMBOL", multiVals = "first")
      )
      keep <- !is.na(fc_map) & !duplicated(fc_map)
      gene_list_fc <- setNames(gene_list_fc[keep], fc_map[keep])
    } else {
      gene_list_fc <- gene_list_fc[!duplicated(names(gene_list_fc))]
    }
  }

  # EnrichR
  enr <- enrichr(symbols, dbs)
  df <- enr[["GO_Biological_Process_2023"]]
  if (is.null(df) || nrow(df) == 0) {
    warning(paste0("run_go_enrichment(", tag, "): Enrichr returned no ",
                    "enriched GO_Biological_Process_2023 terms (empty result, not an error) - enrichR plot skipped."))
  } else {
    df$Count <- as.numeric(gsub("/.*", "", df$Overlap))
    df$Term <- str_wrap(df$Term, width = 30)
    png(file.path(plot_dir, paste0(axis_prefix, "02_GO_", tag, "_enrichR.png")), width = 1000, height = 800, res = 120)
    print(plotEnrich(df, showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value") +
            ggtitle(paste("Top Processes", tag, "-regulated", label)) +
            theme(axis.text.y = element_text(lineheight = 0.8)))
    dev.off()
  }

  # clusterProfiler - BP
  go_bp <- do.call(enrichGO, c(go_args, list(ont = "BP", pvalueCutoff = 0.05)))
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    png(file.path(plot_dir, paste0(axis_prefix, "02_GO_", tag, "_enrichGO.png")), width = 900, height = 700, res = 120)
    print(dotplot(go_bp, showCategory = 10) + ggtitle(paste("Top Processes", tag, "-regulated", label)))
    dev.off()
  } else {
    warning(paste0("run_go_enrichment(", tag, "): enrichGO (BP) found no enriched terms - dotplot skipped."))
  }

  # clusterProfiler - CC
  go_cc <- do.call(enrichGO, c(go_args, list(ont = "CC")))
  if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
    png(file.path(plot_dir, paste0(axis_prefix, "02_GO_", tag, "_CC.png")), width = 1000, height = 800, res = 120)
    print(dotplot(go_cc) + ggtitle(paste("Cellular Components:", tag, "-regulated")))
    dev.off()
  } else {
    warning(paste0("run_go_enrichment(", tag, "): enrichGO (CC) found no enriched terms - dotplot skipped."))
  }

  # clusterProfiler - MF (same design as BP and CC; no significant term on
  # Axis A, reported in the thesis for completeness)
  go_mf <- do.call(enrichGO, c(go_args, list(ont = "MF")))
  if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
    png(file.path(plot_dir, paste0(axis_prefix, "02_GO_", tag, "_MF.png")), width = 1000, height = 800, res = 120)
    print(dotplot(go_mf) + ggtitle(paste("Molecular Function:", tag, "-regulated")))
    dev.off()
  } else {
    warning(paste0("run_go_enrichment(", tag, "): enrichGO (MF) found no enriched terms - dotplot skipped."))
  }

  # Cnetplot
  if (!is.null(gene_list_fc) && !is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    png(file.path(plot_dir, paste0(axis_prefix, "02_GO_", tag, "_Cnetplot.png")), width = 1200, height = 1000, res = 130)
    print(cnetplot(go_bp, foldChange = gene_list_fc, showCategory = 5) +
            ggtitle(paste("Gene-Concept Network:", label)))
    dev.off()
  }

  # Save the complete results as CSV tables, besides the plots, so that the
  # numbers in the thesis come from files. Folder: <output_dir>/tables.
  tables_dir_go <- file.path(dirname(plot_dir), "tables")
  dir.create(tables_dir_go, showWarnings = FALSE, recursive = TRUE)
  save_go_tab <- function(obj, suffix) {
    if (!is.null(obj)) {
      tab <- if (inherits(obj, "enrichResult")) as.data.frame(obj@result) else as.data.frame(obj)
      write.csv(tab, file.path(tables_dir_go, paste0(axis_prefix, "02_GO_", tag, "_", suffix, ".csv")),
                row.names = FALSE)
    }
  }
  save_go_tab(go_bp, "BP_enrichGO")
  save_go_tab(go_cc, "CC_enrichGO")
  save_go_tab(go_mf, "MF_enrichGO")
  if (!is.null(enr)) save_go_tab(enr[["GO_Biological_Process_2023"]], "BP_enrichR")
  if (!is.null(universe_symbols)) {
    writeLines(as.character(length(universe_symbols)),
               file.path(tables_dir_go, paste0(axis_prefix, "02_GO_universe_size.txt")))
  }

  return(list(enrichR = enr, go_bp = go_bp, go_cc = go_cc, go_mf = go_mf))
}


# Functional GSEA (preranked), same analysis for all axes ----
# Complements the ORA of run_go_enrichment(): it uses the ranking of all
# tested genes (Wald statistic for DESeq2, moderated t for limma), so it
# also works with 0 DEGs (Axis B).
# Gene sets: MSigDB Hallmark (50 sets) and GO Biological Process (C5 GO:BP),
# from msigdbr, indexed by Entrez Gene ID. Gene IDs are harmonized with the
# Level 3 dictionary (build_gene_dictionary() + harmonize_deg_table()).
# Outputs (prefix = axis_prefix, e.g. "B"):
#   tables/<P>02_GSEA_Hallmark.csv, tables/<P>02_GSEA_GOBP.csv
#   plots/<P>02_GSEA_Hallmark.png  (NES of all Hallmark sets)
#   plots/<P>02_GSEA_GOBP_top.png  (top GO:BP sets)
#   data/<P>02_GSEA_functional.rds (for the cross-axis comparison)

get_msigdb_sets <- function(collection, subcollection = NULL,
                            cache_dir = "Transferability/Dictionary") {
  # Local copy: after the first successful download the gene sets are saved
  # in cache_dir (msigdb_<collection>.rds), so later runs do not need the
  # network and always use the same MSigDB version.
  tag <- gsub("[^A-Za-z0-9]+", "_", paste(c(collection, subcollection), collapse = "_"))
  cache <- file.path(cache_dir, paste0("msigdb_", tag, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  # msigdbr >= 10 uses collection/subcollection and the ncbi_gene column;
  # older versions use category/subcategory and entrez_gene.
  ms <- tryCatch(
    if (is.null(subcollection)) msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
    else msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection),
    error = function(e) {
      if (is.null(subcollection)) msigdbr::msigdbr(species = "Homo sapiens", category = collection)
      else msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
    })
  id_col <- if ("ncbi_gene" %in% colnames(ms)) "ncbi_gene" else "entrez_gene"
  ms <- ms[!is.na(ms[[id_col]]), ]
  sets <- lapply(split(as.character(ms[[id_col]]), ms$gs_name), unique)
  attr(sets, "msigdbr_version") <- as.character(utils::packageVersion("msigdbr"))
  attr(sets, "downloaded") <- as.character(Sys.Date())
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sets, cache)
  sets
}

# GO Biological Process: MSigDB C5 GO:BP first (needs the msigdbdf data
# package with msigdbr >= 10); otherwise the same sets are built from
# org.Hs.eg.db (propagated GO annotations, GOALL) and GO.db for the term
# names. Set names are "GOBP_<TERM>" in both cases.
get_go_bp_sets <- function() {
  sets <- tryCatch(get_msigdb_sets("C5", "GO:BP"), error = function(e) NULL)
  if (!is.null(sets) && length(sets) > 1000) {
    message("GO:BP from MSigDB (", length(sets), " set)")
    return(sets)
  }
  message("GO:BP from MSigDB not available: using org.Hs.eg.db + GO.db")
  eg <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")
  go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = eg,
                          columns = c("GOALL", "ONTOLOGYALL"), keytype = "ENTREZID"))
  go <- unique(go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", c("ENTREZID", "GOALL")])
  sets <- lapply(split(go$ENTREZID, go$GOALL), unique)
  terms <- AnnotationDbi::Term(GO.db::GOTERM[names(sets)])
  names(sets) <- paste0("GOBP_", toupper(gsub("[^A-Za-z0-9]+", "_", terms[names(sets)])))
  sets[!duplicated(names(sets))]
}

run_functional_gsea <- function(res, axis_prefix, tables_dir, plot_dir, data_dir,
                                stat_col = "stat", dict_dir = "Transferability/Dictionary",
                                min_size = 15, max_size = 500, seed = 18) {
  for (p in c("fgsea", "msigdbr", "GO.db")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("Missing package '", p, "': install.packages(\"msigdbr\") / BiocManager::install(\"fgsea\")")
  }
  dict   <- build_gene_dictionary(dict_dir)
  h      <- harmonize_deg_table(res, dict, stat_col = stat_col)
  ranked <- build_ranked_stats(h, stat_col = stat_col)
  message(axis_prefix, ": ranking of ", length(ranked), " genes (Entrez ID) for the functional GSEA")

  collections <- list(Hallmark = get_msigdb_sets("H"),
                      GOBP     = get_go_bp_sets())
  out <- list()
  for (nm in names(collections)) {
    set.seed(seed)
    r <- fgsea::fgsea(pathways = collections[[nm]], stats = ranked,
                      minSize = min_size, maxSize = max_size, eps = 0,
                      nproc = 1)   # reproducible with set.seed() (see above)
    r <- r[order(r$pval), ]
    r$leadingEdge_symbols <- vapply(r$leadingEdge, function(x)
      paste(na.omit(unname(dict$id2sym[x])), collapse = ";"), character(1))
    r$leadingEdge <- vapply(r$leadingEdge, paste, character(1), collapse = ";")
    r <- as.data.frame(r)
    write.csv(r, file.path(tables_dir, paste0(axis_prefix, "02_GSEA_", nm, ".csv")), row.names = FALSE)
    out[[nm]] <- r
    message(axis_prefix, " - ", nm, ": ", nrow(r), " sets tested, ",
            sum(r$padj < 0.05, na.rm = TRUE), " with padj < 0.05 (",
            sum(r$padj < 0.05 & r$NES > 0, na.rm = TRUE), " NES>0, ",
            sum(r$padj < 0.05 & r$NES < 0, na.rm = TRUE), " NES<0)")
  }

  # Plot 1: NES of all Hallmark sets (filled = padj < 0.05)
  hm <- out$Hallmark
  hm$label <- gsub("_", " ", sub("^HALLMARK_", "", hm$pathway))
  hm$sig <- ifelse(!is.na(hm$padj) & hm$padj < 0.05, "padj < 0.05", "n.s.")
  hm$dir <- ifelse(hm$NES > 0, "Higher in exposed", "Lower in exposed")
  g1 <- ggplot(hm, aes(x = NES, y = reorder(label, NES), fill = dir, alpha = sig)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Higher in exposed" = "#C0392B", "Lower in exposed" = "steelblue"), name = NULL) +
    scale_alpha_manual(values = c("padj < 0.05" = 1, "n.s." = 0.3), name = NULL) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    labs(x = "Normalized enrichment score (NES)", y = NULL) +
    theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
  ggsave(file.path(plot_dir, paste0(axis_prefix, "02_GSEA_Hallmark.png")), g1, width = 6.5, height = 8, dpi = 300)

  # Plot 2: top GO:BP sets (up to 10 per direction, ordered by p)
  gb <- out$GOBP
  top <- rbind(head(gb[gb$NES > 0, ], 10), head(gb[gb$NES < 0, ], 10))
  top$label <- stringr::str_trunc(tolower(gsub("_", " ", sub("^GOBP_", "", top$pathway))), 60)
  g2 <- ggplot(top, aes(x = NES, y = reorder(label, NES), size = size, colour = -log10(padj))) +
    geom_point() + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_colour_gradient(low = "grey70", high = "#8E1B10", name = expression(-log[10]~padj)) +
    labs(x = "Normalized enrichment score (NES)", y = NULL, size = "Set size") +
    theme_minimal(base_size = 9)
  ggsave(file.path(plot_dir, paste0(axis_prefix, "02_GSEA_GOBP_top.png")), g2, width = 7.5, height = 6, dpi = 300)

  out$ranked <- ranked
  saveRDS(out, file.path(data_dir, paste0(axis_prefix, "02_GSEA_functional.rds")))
  invisible(out)
}

# Master pipeline ----

# config: one element of DIMENSIONS (00_config.R). Runs the whole Level 1
# pipeline for that axis and saves everything in <axis folder>/03_results
# (plots, tables, data) and in the SWIMmeR input folder.
run_dimension_pipeline <- function(config) {

  plot_dir <- file.path(config$output_dir, "plots")
  data_dir <- file.path(config$output_dir, "data")
  tables_dir <- file.path(config$output_dir, "tables")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. Import and condition labels
  raw <- load_geo_dataset(config$dataset_id, config$source_dir, config$raw_counts_pattern)
  metadata <- build_condition_variable(raw$metadata,
                                        characteristics_col = config$characteristics_col,
                                        positive_values = config$positive_values)

  ap <- config$axis_prefix  # short prefix for the plot file names

  # 2. Preprocessing (also saves <ap>00_low_expression_filter.png and <ap>00_library_size.png)
  pre <- preprocess_axis_data(raw$raw_counts, metadata, plot_dir = plot_dir, axis_prefix = ap)

  # 3. Batch correction and normalization
  norm <- run_normalization_pipeline(pre$counts, pre$meta,
                                      condition_col = config$condition_col,
                                      group_label = config$group_label,
                                      plot_dir = plot_dir, axis_prefix = ap)

  # 4. DEG identification (also saves <ap>_IQR_distribution.png)
  deg <- run_deg_analysis(norm$dds, norm$vst, condition_col = config$condition_col,
                           positive_level = "Yes", negative_level = "No",
                           plot_dir = plot_dir, axis_prefix = ap)

  # 5. Main plots
  save_volcano_plot(deg$res_filtered, file.path(plot_dir, paste0(ap, "01_volcano.png")), title = "Volcano plot")

  annotation_col <- data.frame(setNames(list(norm$meta[[config$condition_col]]), config$group_label))
  rownames(annotation_col) <- rownames(norm$meta)
  anno_colors <- setNames(list(setNames(c("green", "orange"), c("No", "Yes"))), config$group_label)

  heatmap_matrix <- norm$vst[deg$sig_genes, ]
  save_heatmap_plot(heatmap_matrix, file.path(plot_dir, paste0(ap, "01_heatmap.png")), annotation_col, anno_colors)

  save_boxplot_top_gene(norm$vst, deg$res_filtered, norm$meta, config$condition_col, config$group_label,
                         direction = "up", path = file.path(plot_dir, paste0(ap, "01_boxplot_up.png")))
  save_boxplot_top_gene(norm$vst, deg$res_filtered, norm$meta, config$condition_col, config$group_label,
                         direction = "down", path = file.path(plot_dir, paste0(ap, "01_boxplot_down.png")))

  # 6. Export tables and SWIMmeR input
  export_deg_table(deg$export_tab, file.path(data_dir, "DEG_output.txt"))
  export_top10_xlsx(deg$export_tab, file.path(tables_dir, "Tables_DEGs.xlsx"))
  export_swimmer_files(norm$meta, config$condition_col, counts(norm$dds, normalized = TRUE), config$swimmer_dir)

  # 7. Functional enrichment
  genes_up <- rownames(deg$res_filtered)[deg$res_filtered$diffexpressed == "up"]
  genes_down <- rownames(deg$res_filtered)[deg$res_filtered$diffexpressed == "down"]
  gene_list_fc <- deg$res_filtered$log2FoldChange
  names(gene_list_fc) <- rownames(deg$res_filtered)

  go_down <- run_go_enrichment(genes_down, "down", plot_dir, gene_list_fc, label = config$group_label, axis_prefix = ap)
  go_up <- run_go_enrichment(genes_up, "up", plot_dir, gene_list_fc, label = config$group_label, axis_prefix = ap)

  return(list(metadata = norm$meta, dds = norm$dds, vst = norm$vst,
              deg = deg, go_up = go_up, go_down = go_down))
}

# test_2x2 ----
# Independence test for 2x2 contingency tables: Pearson chi-square without
# Yates correction, or Fisher's exact test when an expected count is below
# 5 (Cochran's rule).
test_2x2 <- function(tb) {
  ct <- suppressWarnings(chisq.test(tb, correct = FALSE))
  if (any(ct$expected < 5)) {
    ft <- fisher.test(tb)
    list(test = "Fisher", statistic = NA, df = NA, p = ft$p.value,
         min_expected = min(ct$expected))
  } else {
    list(test = "chisq", statistic = unname(ct$statistic),
         df = unname(ct$parameter), p = ct$p.value,
         min_expected = min(ct$expected))
  }
}

# hist_pvalues_grid ----
# Histogram of the nominal p-values, with a grey reference grid and a dashed
# line for the uniform distribution expected if no gene differed.
hist_pvalues_grid <- function(p, main, breaks = 50) {
  p <- p[!is.na(p)]
  h <- hist(p, breaks = seq(0, 1, length.out = breaks + 1), plot = FALSE)
  plot(h, col = NA, border = NA, main = main, xlab = "p-value", ylab = "Frequency",
       ylim = c(0, max(h$counts) * 1.25))   # room for the legend
  grid(nx = 10, ny = NULL, col = "grey85", lty = 1)
  plot(h, col = "steelblue", border = "white", add = TRUE)
  abline(h = length(p) / breaks, lty = 2, col = "grey30")
  legend("top", legend = "Uniform distribution (overall mean)", lty = 2, col = "grey30",
         bty = "o", bg = "white", box.col = "grey70", cex = 0.85)
  box(bty = "l")
}
