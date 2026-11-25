################################################################################
#          GENE IDENTIFIER HARMONIZATION (Entrez Gene ID) - ALL AXES
################################################################################
#namefile: Transferability/00_harmonize_gene_ids.R
#creation dt: 08/10/2026
#
# Level 3, step 0. Run AFTER 01_DEG.R of each axis and BEFORE
# Transferability/01_signature_transfer_gsea.R. For each axis it reads the
# complete DEG results table (DEG_full_results.rds; for D the version with gene
# symbols, DEG_full_results_symbols.rds), converts the identifiers to NCBI
# Entrez Gene IDs with build_gene_dictionary() / harmonize_deg_table()
# (00_functions.R) and saves, without touching the axis folders:
#   Transferability/data/<axis>_DEG_full_results_entrez.rds
# plus a summary of the conversion:
#   Transferability/results/harmonization_summary.csv
#
# Why after the DEG analysis: the DEG test of each axis does not depend on gene
# names; names matter only to compare different axes (Level 3). Converting
# earlier would mean merging rows of raw counts.
#
# Run from the project root (open Project.Rproj):
#   source("Transferability/00_harmonize_gene_ids.R", echo = TRUE)

source("00_config.R")
source("00_functions.R")
load_pipeline_libraries()

dict <- build_gene_dictionary("Transferability/Dictionary")
# dictionary sizes
length(dict$sym2id); length(dict$alias2id); length(dict$prev2id); length(dict$ens2id)

inputs <- list(
  A = file.path(DIMENSIONS$mass_trauma$output_dir,     "data", "DEG_full_results.rds"),
  B = file.path(DIMENSIONS$war_trauma$output_dir,      "data", "DEG_full_results.rds"),
  C = file.path(DIMENSIONS$poverty_ses$output_dir,     "data", "DEG_full_results.rds"),
  D = file.path(DIMENSIONS$early_adversity$output_dir, "data", "DEG_full_results_symbols.rds")
)

summary_rows <- list()
harmonized <- list()
for (ax in names(inputs)) {
  f <- inputs[[ax]]
  if (!file.exists(f)) {
    message("Axis ", ax, ": file not found (", f, ") - skipped. Run 01_DEG.R for this axis again.")
    next
  }
  res <- readRDS(f)
  h <- harmonize_deg_table(res, dict)
  out <- file.path("Transferability", "data", paste0(ax, "_DEG_full_results_entrez.rds"))
  saveRDS(h, out)
  harmonized[[ax]] <- h

  s <- attr(h, "mapping_summary")
  summary_rows[[ax]] <- data.frame(
    axis = ax, n_input = nrow(res),
    ensembl = s[["ensembl"]], symbol = s[["symbol"]], alias = s[["alias"]],
    loc_id = s[["loc_id"]], excel_recovered = s[["excel_recovered"]],
    excel_corrupted = s[["excel_corrupted"]],
    ambiguous = s[["ambiguous"]], unmapped = s[["unmapped"]],
    collapsed_duplicates = attr(h, "n_collapsed_duplicates"),
    n_output = nrow(h)
  )
  message("Axis ", ax, ": ", nrow(res), " rows -> ", nrow(h), " genes with an Entrez ID (saved to ", out, ")")
}

summary_df <- do.call(rbind, summary_rows)

# gene names corrupted by Excel, by axis (recovered / excluded)
for (ax in names(harmonized)) {
  h_all <- to_entrez(rownames(readRDS(inputs[[ax]])), dict)
  ex <- h_all[h_all$method %in% c("excel_recovered", "excel_corrupted"), ]
  if (nrow(ex) > 0) {
    ex$symbol_current <- unname(dict$id2sym[ex$entrez])
    message("Axis ", ax, ": gene names corrupted by Excel in the tested table")
    print(ex[order(ex$method, ex$input), ], row.names = FALSE)
    write.csv(ex, file.path("Transferability", "results", paste0("excel_names_", ax, ".csv")), row.names = FALSE)
  }
}
# conversion summary by axis
summary_df
write.csv(summary_df, file.path("Transferability", "results", "harmonization_summary.csv"), row.names = FALSE)

# Objects kept in memory: dict, harmonized (list by axis), summary_df.
