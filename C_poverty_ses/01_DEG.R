################################################################################
#                  DEG - CHRONIC POVERTY (Axis C, GSE221766)
################################################################################
#namefile: 01_DEG.R
#creation dt: 15/08/2026
#
# Axis C, step 2 of 3: DEG identification, main plots and export. Run right
# after 00_preprocessing.R in the same R session: it reuses norm, pre, cfg
# and the output folders already in memory. Functional enrichment is in
# 02_enrichment.R.

# DEG identification (IQR filter + DESeq2 results) ----
# Same thresholds as all axes (FDR 0.05, FC 1.1); the count at FDR < 0.10,
# the threshold of the original study, is printed below. Also saves
# plots/C01_IQR_distribution.png
source("00_functions.R")

deg <- run_deg_analysis(norm$dds, norm$vst, condition_col = cfg$condition_col,
                        positive_level = "Yes", negative_level = "No",
                        plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                        thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc, prc_iqr = cfg$prc_iqr,
                        iqr_plot_xlim = c(0, 1.5))

table(deg$res_filtered$diffexpressed)
length(deg$sig_genes)
deg$export_tab
sum(deg$res_filtered$padj < 0.10, na.rm = TRUE)

# Is XIST (ENSG00000229807) still among the significant genes?
"ENSG00000229807" %in% deg$sig_genes

# Main plots (volcano, p-value histogram, heatmap, boxplots) ----
save_volcano_plot(deg$res_filtered, file.path(plot_dir, paste0(cfg$axis_prefix, "01_volcano.png")), title = "Volcano plot",
                  thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc)

# Histogram of the nominal p-values, as for Axis B
png(file.path(plot_dir, paste0(cfg$axis_prefix, "01_pvalue_histogram.png")),
    width = 800, height = 600, res = 100)
hist_pvalues_grid(deg$res_filtered$pvalue, main = "Distribution of nominal p-values, Axis C")
dev.off()
mean(deg$res_filtered$pvalue < 0.05, na.rm = TRUE)   # share with p < 0.05 (0.05 expected with no effect)
mean(deg$res_filtered$pvalue > 0.9,  na.rm = TRUE)   # share with p > 0.9 (0.10 expected)

annotation_col <- data.frame(setNames(list(norm$meta[[cfg$condition_col]]), cfg$group_label))
rownames(annotation_col) <- rownames(norm$meta)
anno_colors <- setNames(list(setNames(c("green", "orange"), c("No", "Yes"))), cfg$group_label)

heatmap_matrix <- norm$vst[deg$sig_genes, ]
save_heatmap_plot(heatmap_matrix, file.path(plot_dir, paste0(cfg$axis_prefix, "01_heatmap.png")), annotation_col, anno_colors)

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

# Objects kept in memory: raw, metadata, pre, pca, norm, deg (plus cfg and
# the output folders). Next: 02_enrichment.R, in the same R session.
