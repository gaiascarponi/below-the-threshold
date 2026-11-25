################################################################################
#       FUNCTIONAL ENRICHMENT (GO / enrichR) - CHRONIC POVERTY (Axis C)
################################################################################
#namefile: 02_enrichment.R
#creation dt: 15/08/2026
#
# Axis C, step 3 of 3: over-representation analysis (GO, Enrichr) of the
# up- and down-regulated genes and functional GSEA on all tested genes. Run
# right after 00_preprocessing.R and 01_DEG.R in the same R session.

# ORA (GO / enrichR) on up- and down-regulated genes ----
# Plots are saved with the 02_ prefix (plots/C02_GO_...).
genes_up   <- rownames(deg$res_filtered)[deg$res_filtered$diffexpressed == "up"]
genes_down <- rownames(deg$res_filtered)[deg$res_filtered$diffexpressed == "down"]
gene_list_fc <- deg$res_filtered$log2FoldChange
names(gene_list_fc) <- rownames(deg$res_filtered)
# enrichGO background: the genes actually tested (after the IQR filter, non-NA padj), not the whole genome (see run_go_enrichment()).
universe_genes <- rownames(deg$res_filtered)[!is.na(deg$res_filtered$padj)]

go_down <- run_go_enrichment(genes_down, "down", plot_dir, gene_list_fc, label = cfg$group_label, axis_prefix = cfg$axis_prefix,
                              universe = universe_genes)
go_up   <- run_go_enrichment(genes_up, "up", plot_dir, gene_list_fc, label = cfg$group_label, axis_prefix = cfg$axis_prefix,
                              universe = universe_genes)

# Functional GSEA (preranked) on all tested genes ----
# Hallmark and GO:BP, same function for all axes (run_functional_gsea()).
# Ensembl IDs are converted to Entrez inside the function.
gsea_fun <- run_functional_gsea(deg$res_filtered, axis_prefix = cfg$axis_prefix,
                                tables_dir = tables_dir, plot_dir = plot_dir,
                                data_dir = data_dir, stat_col = "stat")
# Saves tables/<P>02_GSEA_*.csv, plots/<P>02_GSEA_*.png and data/<P>02_GSEA_functional.rds
gsea_fun$Hallmark[, c("pathway", "NES", "pval", "padj", "size")]
head(gsea_fun$GOBP[, c("pathway", "NES", "pval", "padj", "size")], 15)

# Objects kept in memory: go_up, go_down, gsea_fun (plus the objects of the
# previous scripts).
