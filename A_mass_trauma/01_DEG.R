################################################################################
#                       DEG - MASS TRAUMA (Axis A, WTC)
################################################################################
#namefile: 01_DEG.R
#creation dt: 18/08/2026
#
# Axis A, step 2 of 3: DEG identification, main plots and export. Run right
# after 00_preprocessing.R in the same R session: it reuses norm, pre, cfg
# and the output folders already in memory. If the session was closed, run
# 00_preprocessing.R again first. Functional enrichment is in 02_enrichment.R.

# DEG identification (IQR filter + DESeq2 results) ----
deg <- run_deg_analysis(norm$dds, norm$vst, condition_col = cfg$condition_col,
                         positive_level = "Yes", negative_level = "No",
                         plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                         thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc, prc_iqr = cfg$prc_iqr,
                        iqr_plot_xlim = c(0, 1.5))
# Saves plots/A01_IQR_distribution.png
table(deg$res_filtered$diffexpressed)
length(deg$sig_genes)
head(deg$export_tab)

# Main plots (volcano, heatmaps, boxplots of the top genes) ----
save_volcano_plot(deg$res_filtered, file.path(plot_dir, paste0(cfg$axis_prefix, "01_volcano.png")), title = "Volcano plot",
                   thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc)
annotation_col <- data.frame(setNames(list(norm$meta[[cfg$condition_col]]), cfg$group_label))
rownames(annotation_col) <- rownames(norm$meta)
anno_colors <- setNames(list(setNames(c("green", "orange"), c("No", "Yes"))), cfg$group_label)

heatmap_matrix <- norm$vst[deg$sig_genes, ]
save_heatmap_plot(heatmap_matrix, file.path(plot_dir, paste0(cfg$axis_prefix, "01_heatmap.png")), annotation_col, anno_colors)

# Second heatmap: samples ordered by group (No, then Yes), without
grp <- norm$meta[[cfg$condition_col]]
ord <- order(factor(grp, levels = c("No", "Yes")))
heatmap_matrix_grouped <- norm$vst[deg$sig_genes, ord]
# column clustering; genes clustered by row.
save_heatmap_plot(heatmap_matrix_grouped, file.path(plot_dir, paste0(cfg$axis_prefix, "01_heatmap_grouped.png")),
                  annotation_col, anno_colors, cluster_cols = FALSE, gaps_col = sum(grp == "No"))

save_boxplot_top_gene(norm$vst, deg$res_filtered, norm$meta, cfg$condition_col, cfg$group_label,
                       direction = "up", path = file.path(plot_dir, paste0(cfg$axis_prefix, "01_boxplot_up.png")))
save_boxplot_top_gene(norm$vst, deg$res_filtered, norm$meta, cfg$condition_col, cfg$group_label,
                       direction = "down", path = file.path(plot_dir, paste0(cfg$axis_prefix, "01_boxplot_down.png")))
# PNG files in plot_dir

# Export: DEG tables, gene lists for Level 2 and SWIMmeR input ----
export_deg_table(deg$export_tab, file.path(data_dir, "DEG_output.txt"))
export_deg_lists(deg$export_tab, file.path(dirname(cfg$output_dir), "04_network", "input"), cfg$axis_prefix)
export_top10_xlsx(deg$export_tab, file.path(tables_dir, "Tables_DEGs.xlsx"))
export_swimmer_files(norm$meta, cfg$condition_col, counts(norm$dds, normalized = TRUE), cfg$swimmer_dir)
# Tables_DEGs.xlsx and DEG_output.txt

# Complete results table (all genes after the IQR filter) --> input of the Level 3 signature transfer
export_deg_full_table(deg$res_filtered, file.path(data_dir, "DEG_full_results.rds"))

# Objects kept in memory: raw, metadata, pre, norm, deg (plus cfg and the
# output folders). Next: 02_enrichment.R, in the same R session.
