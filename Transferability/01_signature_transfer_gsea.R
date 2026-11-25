################################################################################
#      LEVEL 3 - TRANSFER OF THE WTC SIGNATURE (AXIS A) TO B, C, D (GSEA)
################################################################################
#namefile: Transferability/01_signature_transfer_gsea.R
#creation dt: 18/08/2026
#
# The question is not whether B, C and D have significant DEGs of their own
# (they have none), but whether the Axis A signature (525 DEGs, split into up
# and down) is enriched at the top or bottom of the COMPLETE gene ranking of
# each other axis: a GSEA test (fgsea), not an overlap between DEG lists.
#
# Requirements:
#   - DEG_full_results.rds of each axis, saved by 01_DEG.R
#     (export_deg_full_table()); for Axis D, DEG_full_results_symbols.rds,
#     with probe IDs mapped to gene symbols (map_probes_to_symbols_agilent())
#   - Transferability/00_harmonize_gene_ids.R, which writes
#     Transferability/data/<axis>_DEG_full_results_entrez.rds. Genes are
#     matched across axes by NCBI Entrez Gene ID, not by symbol: the datasets
#     use nomenclatures of different years, and matching by symbol lost about
#     70 Axis A genes on Axis C (e.g. CTSL vs CTSL1, KMT2D vs MLL2).
#
# Run from the project root (open Project.Rproj):
#   source("Transferability/01_signature_transfer_gsea.R", echo = TRUE)
# Results in Transferability/results/, figures in Transferability/figures/.

source("00_config.R")
source("00_functions.R")
source("00_source.R")
load_pipeline_libraries()

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea")
}
library(fgsea)

# Source signature (Axis A) and gene sets ----
cfg_A <- DIMENSIONS$mass_trauma
res_A <- readRDS(file.path(cfg_A$output_dir, "data", "DEG_full_results.rds"))
dict  <- build_gene_dictionary("Transferability/Dictionary")

# Signature built from the original DEG list (Axis A symbols), then converted
# to Entrez, so that no DEG is lost when duplicated rows are merged.
up_sym   <- rownames(res_A)[res_A$diffexpressed == "up"]
down_sym <- rownames(res_A)[res_A$diffexpressed == "down"]
map_up   <- to_entrez(up_sym, dict)
map_down <- to_entrez(down_sym, dict)
# how the Axis A DEGs were converted (symbol / alias / loc_id / ...)
table(map_up$method); table(map_down$method)

genes_up_A   <- unique(na.omit(map_up$entrez))
genes_down_A <- unique(na.omit(map_down$entrez))
genes_all_A  <- unique(c(genes_up_A, genes_down_A))

# DEGs converted to Entrez (out of 211 up / 314 down / 525 in total)
length(genes_up_A)
length(genes_down_A)
length(genes_all_A)

signatures_A <- list(
  WTC_up   = genes_up_A,
  WTC_down = genes_down_A,
  WTC_all  = genes_all_A
)

# Axis C (target): ranked vector and test of the Axis A signature ----
cfg_C <- DIMENSIONS$poverty_ses
res_C <- readRDS(file.path("Transferability", "data", "C_DEG_full_results_entrez.rds"))

ranked_C <- build_ranked_stats(res_C)
# genes in the ranked vector and range of the statistic
length(ranked_C)
range(ranked_C)

gsea_C <- run_signature_transfer_gsea(signatures_A, ranked_C)
# pathway, pval, padj, NES (normalized enrichment score) and size (number of
# Axis A signature genes found in the ranked vector)
gsea_C[, c("pathway", "pval", "padj", "NES", "size")]

# Axis D (target) ----
# Uses DEG_full_results_symbols.rds from D_early_adversity/01_DEG.R
# (probe IDs already mapped to gene symbols).
cfg_D <- DIMENSIONS$early_adversity
res_D <- readRDS(file.path("Transferability", "data", "D_DEG_full_results_entrez.rds"))

ranked_D <- build_ranked_stats(res_D)
# genes in the ranked vector and range of the statistic
length(ranked_D)
range(ranked_D)

gsea_D <- run_signature_transfer_gsea(signatures_A, ranked_D)
gsea_D[, c("pathway", "pval", "padj", "NES", "size")]

# Axis B (target, war and combat trauma, GSE64813) ----
cfg_B <- DIMENSIONS$war_trauma
res_B <- readRDS(file.path("Transferability", "data", "B_DEG_full_results_entrez.rds"))
ranked_B <- build_ranked_stats(res_B)
# genes in the ranked vector and range of the statistic
length(ranked_B)
range(ranked_B)
gsea_B <- run_signature_transfer_gsea(signatures_A, ranked_B)
gsea_B[, c("pathway", "pval", "padj", "NES", "size")]

# Save the GSEA results (leading edges translated to symbols) ----
le_to_symbols <- function(g) {
  g$leadingEdge_symbols <- vapply(g$leadingEdge, function(x)
    paste(unname(dict$id2sym[x]), collapse = ";"), character(1))
  g$leadingEdge <- vapply(g$leadingEdge, paste, character(1), collapse = ";")
  g
}
gsea_all <- list(A_to_C = le_to_symbols(gsea_C), A_to_D = le_to_symbols(gsea_D),
                 A_to_B = le_to_symbols(gsea_B))
res_dir_L3 <- file.path("Transferability", "results")
dir.create(res_dir_L3, showWarnings = FALSE, recursive = TRUE)
saveRDS(gsea_all, file.path(res_dir_L3, "gsea_signature_transfer.rds"))
for (nm in names(gsea_all)) {
  write.csv(gsea_all[[nm]], file.path(res_dir_L3, paste0("gsea_", nm, ".csv")), row.names = FALSE)
}

# Enrichment plots (one figure per signature and target axis) ----
plot_dir_L3 <- file.path("Transferability", "figures")
dir.create(plot_dir_L3, showWarnings = FALSE, recursive = TRUE)

targets <- list(B = list(ranked = ranked_B, label = "Axis B (war/combat trauma)"),
                C = list(ranked = ranked_C, label = "Axis C (poverty/SES)"),
                D = list(ranked = ranked_D, label = "Axis D (early adversity)"))
sig_labels <- c(WTC_down = "down-regulated", WTC_up = "up-regulated", WTC_all = "complete")

for (tg in names(targets)) {
  for (sg in names(signatures_A)) {
    p <- plotEnrichment(signatures_A[[sg]], targets[[tg]]$ranked) +
      ggtitle(paste0("WTC ", sig_labels[[sg]], " signature (Axis A) on ", targets[[tg]]$label))
    ggsave(file.path(plot_dir_L3, paste0("L3_Ato", tg, "_", sub("_", "", sg), ".png")),
           p, width = 7, height = 4.5, dpi = 300)
  }
}
# 9 PNG files in Transferability/figures/ (L3_AtoB_WTCdown.png, ...)
list.files(plot_dir_L3, pattern = "^L3_")

# Sensitivity: Axis A signature without histone and hemoglobin genes ----
# The leading edges of A->B and A->D contain many histone and hemoglobin
# genes, which are sensitive to cell composition (red cells, reticulocytes,
# proliferating cells) more than to gene regulation. The list is defined A
# PRIORI from the official HGNC groups (not chosen by looking at the
# results): "Hemoglobin subunits" and "H1/H2A/H2B/H3/H4 histones".
hgnc_tab <- read.delim(file.path("Transferability", "Dictionary", "hgnc_complete_set.txt"),
                       quote = "", check.names = FALSE, colClasses = "character")
hgnc_tab$gene_group <- gsub('"', "", hgnc_tab$gene_group, fixed = TRUE)
excl_groups <- c("Hemoglobin subunits", "H1 histones", "H2A histones",
                 "H2B histones", "H3 histones", "H4 histones")
in_excl <- vapply(strsplit(hgnc_tab$gene_group, "|", fixed = TRUE),
                  function(g) any(g %in% excl_groups), logical(1))
excl_entrez <- unique(hgnc_tab$entrez_id[in_excl & hgnc_tab$entrez_id != ""])

signatures_A_noHH <- lapply(signatures_A, function(g) setdiff(g, excl_entrez))
names(signatures_A_noHH) <- paste0(names(signatures_A), "_noHistHb")
# genes removed from each signature (about 36 expected, all from the down part)
sapply(signatures_A, length) - sapply(signatures_A_noHH, length)

gsea_sens <- list(
  A_to_B = run_signature_transfer_gsea(signatures_A_noHH, ranked_B),
  A_to_C = run_signature_transfer_gsea(signatures_A_noHH, ranked_C),
  A_to_D = run_signature_transfer_gsea(signatures_A_noHH, ranked_D)
)
gsea_sens <- lapply(gsea_sens, le_to_symbols)
# comparison with the main results above
lapply(gsea_sens, function(g) g[, c("pathway", "pval", "padj", "NES", "size")])

saveRDS(gsea_sens, file.path(res_dir_L3, "gsea_sensitivity_noHistHb.rds"))
for (nm in names(gsea_sens)) {
  write.csv(gsea_sens[[nm]], file.path(res_dir_L3, paste0("gsea_sensitivity_noHistHb_", nm, ".csv")),
            row.names = FALSE)
}

# Objects kept in memory: res_A, signatures_A, ranked_B/C/D, gsea_B/C/D,
# gsea_all, dict.
