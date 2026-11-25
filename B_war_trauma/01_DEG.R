################################################################################
#                  DEG - WAR/COMBAT TRAUMA (Axis B, GSE64813)
################################################################################
#namefile: 01_DEG.R
#creation dt: 18/08/2026
#last edit dt: 8/10/2026
#
# Axis B, step 2 of 3: DEG identification, p-value diagnostics, main plots
# and export. Run right after 00_preprocessing.R in the same R session.
# No gene passes the threshold on this axis, so the p-value distribution is
# checked below. The original study (Breen et al., 2015) described an innate
# immune and interferon module, compared with the functional GSEA in
# 02_enrichment.R.

# Imports & Libraries ----
source("00_config.R")
source("00_functions.R")
source("00_source.R")

# DEG identification (IQR filter + DESeq2 results) ----
deg <- run_deg_analysis(norm$dds, norm$vst, condition_col = cfg$condition_col,
                         positive_level = "Yes", negative_level = "No",
                         plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                         thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc, prc_iqr = cfg$prc_iqr,
                        iqr_plot_xlim = c(0, 1.5))
table(deg$res_filtered$diffexpressed)
length(deg$sig_genes)
head(deg$export_tab)

r <- deg$res_filtered
sum(!is.na(r$padj))                       # genes with a padj (after independent filtering)
min(r$padj, na.rm = TRUE)                 # smallest padj
sum(r$pvalue < 0.05, na.rm = TRUE)        # nominal p < 0.05 (about 5% expected by chance)
head(r[order(r$pvalue), c("log2FoldChange","pvalue","padj")], 10)   # top 10 genes
## P-value diagnostics ----
# Histogram of the nominal p-values, saved for the thesis
png(file.path(plot_dir, paste0(cfg$axis_prefix, "01_pvalue_histogram.png")),
    width = 800, height = 600, res = 100)
hist_pvalues_grid(r$pvalue, main = "Distribution of nominal p-values, Axis B")
dev.off()

# Does the conservative shape depend on lowly expressed genes?
par(mfrow = c(1, 2))
hist(r$pvalue[r$baseMean <  50], breaks = 50, main = "baseMean < 50",  xlab = "p-value")
hist(r$pvalue[r$baseMean >= 50], breaks = 50, main = "baseMean >= 50", xlab = "p-value")
par(mfrow = c(1, 1))
table(r$baseMean >= 50)
tapply(r$pvalue < 0.05, r$baseMean >= 50, mean)
# Robustness to the baseMean threshold (thesis, Appendix A): the conservative
# shape should hold for well-expressed genes at any threshold
for (thr in c(20, 50, 100, 200, median(r$baseMean))) {
  hi <- r$baseMean >= thr
  cat(sprintf("threshold %6.1f: n low=%5d high=%5d | p<0.05 low %.1f%% high %.1f%% | p>0.9 low %.1f%% high %.1f%%\n",
              thr, sum(!hi), sum(hi),
              100 * mean(r$pvalue[!hi] < 0.05, na.rm = TRUE), 100 * mean(r$pvalue[hi] < 0.05, na.rm = TRUE),
              100 * mean(r$pvalue[!hi] > 0.9,  na.rm = TRUE), 100 * mean(r$pvalue[hi] > 0.9,  na.rm = TRUE)))
}

png(file.path(plot_dir, paste0(cfg$axis_prefix, "01_pvalue_by_basemean.png")),
    width = 1000, height = 500, res = 100)
par(mfrow = c(1, 2))
hist_pvalues_grid(r$pvalue[r$baseMean <  50], main = "baseMean < 50")
hist_pvalues_grid(r$pvalue[r$baseMean >= 50], main = "baseMean >= 50")
par(mfrow = c(1, 1))
dev.off()

table(r$baseMean >= 50)

# Main plots (volcano, heatmap, boxplots of the top genes) ----
save_volcano_plot(deg$res_filtered, file.path(plot_dir, paste0(cfg$axis_prefix, "01_volcano.png")), title = "Volcano plot",
                   thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc)

annotation_col <- data.frame(setNames(list(norm$meta[[cfg$condition_col]]), cfg$group_label))
rownames(annotation_col) <- rownames(norm$meta)
anno_colors <- setNames(list(setNames(c("green", "orange"), c("No", "Yes"))), cfg$group_label)

heatmap_matrix <- norm$vst[deg$sig_genes, ]
if (length(deg$sig_genes) == 0) {
  message("No gene passes the threshold (empty sig_genes) - heatmap skipped.")
} else {
  save_heatmap_plot(heatmap_matrix, file.path(plot_dir, paste0(cfg$axis_prefix, "01_heatmap.png")), annotation_col, anno_colors)
}

save_boxplot_top_gene(norm$vst, deg$res_filtered, norm$meta, cfg$condition_col, cfg$group_label,
                       direction = "up", path = file.path(plot_dir, paste0(cfg$axis_prefix, "01_boxplot_up.png")))
save_boxplot_top_gene(norm$vst, deg$res_filtered, norm$meta, cfg$condition_col, cfg$group_label,
                       direction = "down", path = file.path(plot_dir, paste0(cfg$axis_prefix, "01_boxplot_down.png")))

# Export: DEG tables, gene lists for Level 2 and SWIMmeR input ----
export_deg_table(deg$export_tab, file.path(data_dir, "DEG_output.txt"))
export_deg_lists(deg$export_tab, file.path(dirname(cfg$output_dir), "04_network", "input"), cfg$axis_prefix)
export_top10_xlsx(deg$export_tab, file.path(tables_dir, "Tables_DEGs.xlsx"))
export_swimmer_files(norm$meta, cfg$condition_col, counts(norm$dds, normalized = TRUE), cfg$swimmer_dir)

# Complete results table (all genes after the IQR filter) --> input of the Level 3 signature transfer
export_deg_full_table(deg$res_filtered, file.path(data_dir, "DEG_full_results.rds"))

# Objects kept in memory: raw, metadata, pre, norm, deg (plus cfg and the
# output folders). Next: 02_enrichment.R, in the same R session.
