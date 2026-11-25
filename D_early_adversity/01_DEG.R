################################################################################
#                DEG - EARLY-LIFE ADVERSITY (Axis D, GSE70603)
################################################################################
#namefile: 01_DEG.R
#creation dt: 15/08/2026
#
# Axis D, step 2 of 3: DEG identification (limma), main plots and export.
# Run right after 00_preprocessing.R in the same R session: it reuses norm,
# pre, cfg and the output folders already in memory. Probe IDs are mapped
# to gene symbols here (res_D_symbols), for Level 3 and 02_enrichment.R.

# DEG identification (limma: IQR filter + topTable) ----
source("00_functions.R")
# Sex as a covariate: with the PCA computed as for the RNA-seq axes (500
# most variable probes, not scaled) sex explains about 99% of the PC1
# variance, as on Axis C; the original authors also adjusted for sex. Sex is
# balanced between groups (20 F / 10 M in each).
norm$meta$sex <- factor(norm$meta[["gender:ch1"]])
table(norm$meta$sex, norm$meta[[cfg$condition_col]])   # 20/10 expected in each group
deg <- run_deg_analysis_microarray(norm$exprs, norm$meta, condition_col = cfg$condition_col,
                                    positive_level = "Yes", negative_level = "No",
                                    plot_dir = plot_dir, axis_prefix = cfg$axis_prefix,
                                    thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc, prc_iqr = cfg$prc_iqr,
                                   iqr_plot_xlim = c(0, 2.5), 
                                   extra_covariates = "sex")
length(deg$iqr_values); sum(deg$iqr_values > deg$thr_perc)   # probes in total / after the IQR filter
min(deg$res_filtered$padj, na.rm = TRUE)                      # smallest padj
table(deg$res_filtered$diffexpressed)
length(deg$sig_genes)

# Main plots (volcano, p-value histogram, boxplots) ----
save_volcano_plot(deg$res_filtered, file.path(plot_dir, paste0(cfg$axis_prefix, "01_volcano.png")), title = "Volcano plot",
                   thr_pval = cfg$thr_pval, thr_fc = cfg$thr_fc)

# Histogram of the nominal p-values, as for Axis B
png(file.path(plot_dir, paste0(cfg$axis_prefix, "01_pvalue_histogram.png")),
    width = 800, height = 600, res = 100)
hist_pvalues_grid(deg$res_filtered$pvalue, main = "Distribution of nominal p-values, Axis D")
dev.off()
mean(deg$res_filtered$pvalue < 0.05, na.rm = TRUE)   # share with p < 0.05 (0.05 expected with no effect)
mean(deg$res_filtered$pvalue > 0.9,  na.rm = TRUE)   # share with p > 0.9 (0.10 expected)

annotation_col <- data.frame(setNames(list(norm$meta[[cfg$condition_col]]), cfg$group_label))
rownames(annotation_col) <- rownames(norm$meta)
anno_colors <- setNames(list(setNames(c("green", "orange"), c("No", "Yes"))), cfg$group_label)

heatmap_matrix <- norm$exprs[deg$sig_genes, ]
#save_heatmap_plot(heatmap_matrix, file.path(plot_dir, paste0(cfg$axis_prefix, "01_heatmap.png")), annotation_col, anno_colors)

save_boxplot_top_gene(norm$exprs, deg$res_filtered, norm$meta, cfg$condition_col, cfg$group_label,
                       direction = "up", path = file.path(plot_dir, paste0(cfg$axis_prefix, "01_boxplot_up.png")))
save_boxplot_top_gene(norm$exprs, deg$res_filtered, norm$meta, cfg$condition_col, cfg$group_label,
                       direction = "down", path = file.path(plot_dir, paste0(cfg$axis_prefix, "01_boxplot_down.png")))

table(deg$res_filtered$diffexpressed)
summary(deg$res_filtered$padj)
range(deg$res_filtered$log2FoldChange)

## Check: same contrast without the sex covariate ----
fit_check <- lmFit(norm$exprs, model.matrix(~ norm$meta[[cfg$condition_col]]))
fit_check <- eBayes(fit_check)
tt_check <- topTable(fit_check, coef = 2, number = Inf)
summary(tt_check$P.Value)
sum(tt_check$P.Value < 0.05, na.rm = TRUE)

# Export: DEG tables, gene lists for Level 2 and SWIMmeR input ----
export_deg_table(deg$export_tab, file.path(data_dir, "DEG_output.txt"))
export_deg_lists(deg$export_tab, file.path(dirname(cfg$output_dir), "04_network", "input"), cfg$axis_prefix)
export_top10_xlsx(deg$export_tab, file.path(tables_dir, "Tables_DEGs.xlsx"))

# Complete results table (all probes after the IQR filter; rownames are
# still Agilent probe IDs)
export_deg_full_table(deg$res_filtered, file.path(data_dir, "DEG_full_results.rds"))

# Same table with probe IDs mapped to gene symbols (GPL17077 annotation):
# the one used for Level 3 (build_ranked_stats()) and in 02_enrichment.R
res_D_symbols <- map_probes_to_symbols_agilent(deg$res_filtered, platform_id = "GPL17077")
# share of mapped probes (printed by the function) and the resulting table
dim(res_D_symbols)
head(res_D_symbols)
export_deg_full_table(res_D_symbols, file.path(data_dir, "DEG_full_results_symbols.rds"))

# SWIMmeR needs a linear-scale matrix; norm$exprs is log2 (limma / ComBat),
# so it is back-transformed with 2^x before the export.
export_swimmer_files(norm$meta, cfg$condition_col, 2^(norm$exprs), cfg$swimmer_dir)


head(res_D_symbols[order(res_D_symbols$pvalue), ], 12)

# Objects kept in memory: raw, metadata, pre, norm, deg, res_D_symbols (plus
# cfg and the output folders). Next: 02_enrichment.R, in the same R session.
