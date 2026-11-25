################################################################################
#            CELL COMPOSITION (DECONVOLUTION) AND THE A->C TRANSFER
################################################################################
#namefile: Transferability/04_cell_composition.R
#creation dt: 14/10/2026
#
# Does the transfer of the Axis A signature to Axis C reflect adversity or
# differences in the cell composition of the samples? (The A->C leading edge
# contains myeloid and erythroid markers: FCN1, AIF1, S100A6, EPOR, ICAM4,
# BLVRB; see 03_leading_edge.R.)
#
# Method
#   1. Cell-type abundance scores with the MCP-counter method (Becht et al.
#      2016, Genome Biol 17:218): score = mean log-normalized expression of the
#      marker genes of each population. Markers from the package file (local
#      copy: Transferability/Dictionary/MCPcounter_genes.txt, MCPcounter
#      v1.2.0, GPL-3), blood populations only (endothelium and fibroblasts
#      excluded). Two extra scores with the same logic, not in MCP-counter:
#      erythroid (hemoglobins and red-cell genes) and platelet.
#      Scores are compared across samples, not across populations.
#   2. Scores compared between the groups of each axis (A, B, C): Wilcoxon
#      and linear model (C: adjusted for sex).
#   3. Axis C: DESeq2 rerun with the cell scores as covariates (same samples,
#      same counts, same IQR filter), new ranking, new A->C transfer. Two
#      pre-specified models:
#        reduced = ~ sex + monocytic + erythroid + condition_bin
#                  (the two components suggested by the leading edge)
#        full    = ~ sex + T + B + NK + monocytic + neutrophils + erythroid +
#                  platelets + condition_bin
#   4. Axis D: NO deconvolution. The samples are CD14+ monocytes isolated with
#      MACS (GEO GSE70603), so the composition is fixed by design: the
#      confounding that deconvolution corrects does not arise, and marker-based
#      scores would only measure residual contamination. A PURITY CHECK is done
#      instead: erythroid, platelet, T, B, NK and monocytic scores compared
#      between groups.
#
# Scores are computed on the VST matrix (RNA-seq) or the quantile-normalized
# matrix (microarray) produced by the 00_preprocessing.R script of each axis,
# which is rerun here in a separate environment. The preprocessing figures are
# regenerated unchanged.
#
# Run from the project root, AFTER 00_harmonize_gene_ids.R and
# 01_signature_transfer_gsea.R:
#   source("Transferability/04_cell_composition.R", echo = TRUE)
# Time: a few minutes (preprocessing of 4 axes + 2 DESeq2 runs on C).
# Outputs: Transferability/results/cell_*.csv, Transferability/figures/cell_*.png/.pdf
################################################################################

source("00_config.R")
source("00_functions.R")
load_pipeline_libraries()
library(fgsea)
set.seed(18)

res_dir <- "Transferability/results"
fig_dir <- "Transferability/figures"

# Markers ----
mcp <- read.delim("Transferability/Dictionary/MCPcounter_genes.txt", check.names = FALSE,
                  stringsAsFactors = FALSE, colClasses = "character")
names(mcp) <- c("symbol", "population", "entrez", "ensembl")
blood_pops <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage",
                "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils")
mcp <- mcp[mcp$population %in% blood_pops, ]
# extra scores (same logic: mean of the markers)
extra <- rbind(
  data.frame(symbol = c("HBB", "HBA1", "HBA2", "HBD", "ALAS2", "CA1", "SLC4A1", "AHSP"),
             population = "Erythroid"),
  data.frame(symbol = c("PPBP", "PF4", "ITGA2B", "GP9", "TUBB1"), population = "Platelets"))
markers_sym <- c(split(mcp$symbol, mcp$population), split(extra$symbol, extra$population))
sapply(markers_sym, length)

# matrix by symbol: rows = symbols (mean of duplicates), columns = samples
to_symbol_matrix <- function(mat, id2sym) {
  sym <- unname(id2sym[rownames(mat)])
  keep <- !is.na(sym) & sym != ""
  rs  <- rowsum(mat[keep, , drop = FALSE], sym[keep])
  cnt <- table(sym[keep])
  rs / as.vector(cnt[rownames(rs)])
}
id2symbol_auto <- function(ids) {
  if (mean(grepl("^ENSG", ids)) > 0.5) {
    core <- sub("\\..*$", "", ids)
    s <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(core), keytype = "ENSEMBL",
                                                column = "SYMBOL", multiVals = "first"))
    setNames(unname(s[core]), ids)
  } else setNames(ids, ids)
}
cell_scores <- function(sym_mat) {
  sc <- sapply(markers_sym, function(g) {
    g <- intersect(g, rownames(sym_mat))
    if (length(g) == 0) return(rep(NA_real_, ncol(sym_mat)))
    colMeans(sym_mat[g, , drop = FALSE])
  })
  sc <- matrix(sc, nrow = ncol(sym_mat), dimnames = list(colnames(sym_mat), names(markers_sym)))
  attr(sc, "n_markers_found") <- sapply(markers_sym, function(g) length(intersect(g, rownames(sym_mat))))
  sc
}

# Preprocessing of each axis (separate environment) and scores ----
run_preproc <- function(script) {
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env, chdir = FALSE)
  env
}
axes <- list(A = "A_mass_trauma/00_preprocessing.R", B = "B_war_trauma/00_preprocessing.R",
             C = "C_poverty_ses/00_preprocessing.R", D = "D_early_adversity/00_preprocessing.R")
env <- list(); scores <- list(); meta <- list()
for (ax in names(axes)) {
  cat("\n================ Axis", ax, "- preprocessing ================\n")
  env[[ax]] <- run_preproc(axes[[ax]])
  n <- env[[ax]]$norm
  if (ax == "D") {
    gpl <- GEOquery::Table(env[[ax]]$gpl_D)
    sc <- intersect(c("Gene symbol", "GENE_SYMBOL", "Symbol"), colnames(gpl))[1]
    id2sym <- setNames(trimws(sub("///.*$", "", as.character(gpl[[sc]]))), as.character(gpl$ID))
    mat <- n$exprs
  } else {
    id2sym <- id2symbol_auto(rownames(n$vst)); mat <- n$vst
  }
  sm <- to_symbol_matrix(as.matrix(mat), id2sym)
  scores[[ax]] <- cell_scores(sm)
  meta[[ax]]   <- n$meta
  cat("Markers found per population (axis", ax, "):\n"); print(attr(scores[[ax]], "n_markers_found"))
}

# Scores compared between groups ----
cond_col <- function(ax) DIMENSIONS[[c(A = "mass_trauma", B = "war_trauma", C = "poverty_ses",
                                       D = "early_adversity")[ax]]]$condition_col
test_scores <- do.call(rbind, lapply(names(scores), function(ax) {
  S <- scores[[ax]]; M <- meta[[ax]][rownames(S), , drop = FALSE]
  grp <- factor(M[[cond_col(ax)]], levels = c("No", "Yes"))
  if (nlevels(droplevels(grp)) < 2) {
    warning("Axis ", ax, ": groups not found as 'No'/'Yes' in ", cond_col(ax), " - check")
    return(NULL)
  }
  do.call(rbind, lapply(colnames(S), function(p) {
    x <- S[, p]; if (all(is.na(x))) return(NULL)
    w <- suppressWarnings(wilcox.test(x ~ grp))
    f <- if (ax == "C" && "sex" %in% names(M)) lm(x ~ M$sex + grp) else lm(x ~ grp)
    cf <- summary(f)$coefficients; r <- grep("^grp", rownames(cf))
    data.frame(axis = ax, population = p, n_markers = attr(S, "n_markers_found")[p],
               mean_No = mean(x[grp == "No"]), mean_Yes = mean(x[grp == "Yes"]),
               diff_Yes_minus_No = cf[r, 1], lm_p = cf[r, 4], wilcox_p = w$p.value)
  }))
}))
test_scores$lm_padj <- ave(test_scores$lm_p, test_scores$axis, FUN = function(p) p.adjust(p, "BH"))
print(test_scores, digits = 3, row.names = FALSE)
write.csv(test_scores, file.path(res_dir, "cell_scores_group_tests.csv"), row.names = FALSE)
for (ax in names(scores))
  write.csv(data.frame(sample = rownames(scores[[ax]]), scores[[ax]], check.names = FALSE),
            file.path(res_dir, paste0("cell_scores_", ax, ".csv")), row.names = FALSE)

# figure: standardized scores by group, one panel per population
long <- do.call(rbind, lapply(names(scores), function(ax) {
  S <- scale(scores[[ax]]); M <- meta[[ax]][rownames(S), , drop = FALSE]
  do.call(rbind, lapply(colnames(S), function(p)
    data.frame(axis = ax, population = p, z = S[, p],
               group = factor(M[[cond_col(ax)]], levels = c("No", "Yes")))))
}))
long <- long[!is.na(long$z), ]
g <- ggplot(long, aes(axis, z, fill = group)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3, width = 0.7) +
  facet_wrap(~ population, nrow = 2) +
  scale_fill_manual(values = c(No = "#9ec5f4", Yes = "#2a78d6"), name = "Exposed") +
  labs(x = "Axis", y = "Cell-type score (z, within axis)",
       title = "Estimated cell-type abundance by group",
       caption = "MCP-counter marker scores (Becht et al. 2016) plus erythroid and platelet scores. Axis D: purified CD14+ monocytes (purity check only).") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "cell_scores_by_group.png"), g, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "cell_scores_by_group.pdf"), g, width = 10, height = 5.5)

# Axis C: DESeq2 with cell covariates and new A->C transfer ----
cfgC  <- DIMENSIONS$poverty_ses
normC <- env$C$norm
SC <- scores$C[colnames(normC$dds), , drop = FALSE]
cd <- as.data.frame(colData(normC$dds))
cd$sizeFactor <- NULL                       # re-estimated by DESeq() (same counts -> same values)
zc <- function(x) as.numeric(scale(x))
cd$cs_T   <- zc(SC[, "T cells"]);     cd$cs_B    <- zc(SC[, "B lineage"])
cd$cs_NK  <- zc(SC[, "NK cells"]);    cd$cs_mono <- zc(SC[, "Monocytic lineage"])
cd$cs_neu <- zc(SC[, "Neutrophils"]); cd$cs_ery  <- zc(SC[, "Erythroid"])
cd$cs_plt <- zc(SC[, "Platelets"])
models <- list(
  reduced = c("sex", "cs_mono", "cs_ery"),
  full    = c("sex", "cs_T", "cs_B", "cs_NK", "cs_mono", "cs_neu", "cs_ery", "cs_plt"))

# Axis A signature (as in 01_signature_transfer_gsea.R)
dict  <- build_gene_dictionary("Transferability/Dictionary")
res_A <- readRDS(file.path(DIMENSIONS$mass_trauma$output_dir, "data", "DEG_full_results.rds"))
sig_up   <- unique(na.omit(to_entrez(rownames(res_A)[res_A$diffexpressed == "up"],   dict)$entrez))
sig_down <- unique(na.omit(to_entrez(rownames(res_A)[res_A$diffexpressed == "down"], dict)$entrez))
signatures_A <- list(WTC_up = sig_up, WTC_down = sig_down, WTC_all = unique(c(sig_up, sig_down)))

ranked_orig <- build_ranked_stats(readRDS("Transferability/data/C_DEG_full_results_entrez.rds"))
gsea_rows <- list(cbind(model = "original (~ sex + condition)",
                        run_signature_transfer_gsea(signatures_A, ranked_orig)[, c("pathway", "NES", "pval", "padj", "size")]))
stat_cor <- c()
for (m in names(models)) {
  cat("\n---- Axis C, model", m, ":", paste(c(models[[m]], cfgC$condition_col), collapse = " + "), "\n")
  dds_m <- DESeqDataSetFromMatrix(countData = counts(normC$dds), colData = cd,
                                  design = as.formula(paste("~", paste(c(models[[m]], cfgC$condition_col), collapse = " + "))))
  dds_m <- DESeq(dds_m, quiet = TRUE)
  deg_m <- run_deg_analysis(dds_m, normC$vst, condition_col = cfgC$condition_col,
                            positive_level = "Yes", negative_level = "No",
                            thr_pval = cfgC$thr_pval, thr_fc = cfgC$thr_fc, prc_iqr = cfgC$prc_iqr)
  cat("DEGs with model", m, ":", length(deg_m$sig_genes), "\n")
  res_m <- harmonize_deg_table(deg_m$res_filtered, dict)
  saveRDS(res_m, file.path("Transferability", "data", paste0("C_DEG_full_results_entrez_celladj_", m, ".rds")))
  ranked_m <- build_ranked_stats(res_m)
  common <- intersect(names(ranked_m), names(ranked_orig))
  stat_cor[m] <- cor(ranked_m[common], ranked_orig[common], method = "spearman")
  g_m <- run_signature_transfer_gsea(signatures_A, ranked_m)
  gsea_rows[[m]] <- cbind(model = paste0(m, " (", paste(models[[m]], collapse = " + "), " + condition)"),
                          g_m[, c("pathway", "NES", "pval", "padj", "size")])
}
gsea_cmp <- do.call(rbind, gsea_rows)
cat("\nSpearman correlation between the original and the adjusted statistic:\n"); print(round(stat_cor, 3))
print(gsea_cmp, digits = 3, row.names = FALSE)
write.csv(gsea_cmp, file.path(res_dir, "cell_adjusted_A_to_C_gsea.csv"), row.names = FALSE)

# Axis D: purity check (summary) ----
cat("\nAxis D (isolated CD14+ monocytes): group differences in contamination scores\n")
print(test_scores[test_scores$axis == "D", c("population", "n_markers", "mean_No", "mean_Yes", "lm_p", "lm_padj")],
      digits = 3, row.names = FALSE)
