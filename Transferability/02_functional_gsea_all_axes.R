################################################################################
#            FUNCTIONAL GSEA (Hallmark + GO:BP) - ALL AXES AT ONCE
################################################################################
#namefile: Transferability/02_functional_gsea_all_axes.R
#creation dt: 13/10/2026
#
# Repeats the functional GSEA of 02_enrichment.R (run_functional_gsea(),
# 00_functions.R) for A, B, C and D, reading the complete DEG results tables
# saved by 01_DEG.R, without rerunning preprocessing and DEG analysis. It
# regenerates all functional GSEAs with the same code (fgsea with nproc = 1,
# reproducible) and the same local copy of MSigDB. Outputs go to each axis
# folder, with the same names as in 02_enrichment.R (tables/<P>02_GSEA_*.csv,
# plots/<P>02_GSEA_*.png, data/<P>02_GSEA_functional.rds), plus a summary in
# Transferability/results/functional_gsea_summary.csv
#
# Run from the project root, AFTER 01_DEG.R of all axes:
#   source("Transferability/02_functional_gsea_all_axes.R", echo = TRUE)

source("00_config.R")
source("00_functions.R")
load_pipeline_libraries()

axes <- list(
  A = list(cfg = DIMENSIONS$mass_trauma,     file = "DEG_full_results.rds"),
  B = list(cfg = DIMENSIONS$war_trauma,      file = "DEG_full_results.rds"),
  C = list(cfg = DIMENSIONS$poverty_ses,     file = "DEG_full_results.rds"),
  D = list(cfg = DIMENSIONS$early_adversity, file = "DEG_full_results_symbols.rds")
)

gsea_all <- list(); summary_rows <- list()
for (ax in names(axes)) {
  cfg_ax <- axes[[ax]]$cfg
  out_dir <- cfg_ax$output_dir
  f <- file.path(out_dir, "data", axes[[ax]]$file)
  if (!file.exists(f)) { message("Axis ", ax, ": missing ", f, " - skipped"); next }
  res <- readRDS(f)
  # file read, last modified date and number of rows
  # (expected: A 19637, B 16764, C 28056, D 26639; see harmonization_summary.csv)
  message("Axis ", ax, ": ", f, " | modified ", format(file.mtime(f), "%d/%m/%Y %H:%M"),
          " | ", nrow(res), " genes")
  g <- run_functional_gsea(res, axis_prefix = cfg_ax$axis_prefix,
                           tables_dir = file.path(out_dir, "tables"),
                           plot_dir   = file.path(out_dir, "plots"),
                           data_dir   = file.path(out_dir, "data"),
                           stat_col   = "stat")
  gsea_all[[ax]] <- g
  for (coll in c("Hallmark", "GOBP")) {
    r <- g[[coll]]; sig <- !is.na(r$padj) & r$padj < 0.05
    summary_rows[[paste(ax, coll)]] <- data.frame(axis = ax, collection = coll,
      n_ranked = length(g$ranked), n_sets = nrow(r), n_sig = sum(sig),
      n_sig_pos = sum(sig & r$NES > 0), n_sig_neg = sum(sig & r$NES < 0))
  }
}
summary_df <- do.call(rbind, summary_rows)
# summary by axis and collection
summary_df
write.csv(summary_df, file.path("Transferability", "results", "functional_gsea_summary.csv"), row.names = FALSE)

# Objects kept in memory: gsea_all (list by axis), summary_df.
