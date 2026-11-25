################################################################################
#      FUNCTIONAL ENRICHMENT (ORA + GSEA) - EARLY-LIFE ADVERSITY (Axis D)
################################################################################
#namefile: 02_enrichment.R
#creation dt: 15/08/2026
#
# Axis D, step 3 of 3: over-representation analysis (GO, Enrichr) and
# functional GSEA. Run right after 00_preprocessing.R and 01_DEG.R in the
# same R session. Both analyses use res_D_symbols (01_DEG.R), where probe
# IDs are mapped to gene symbols (one probe per symbol, largest |t|).

# ORA (GO / enrichR) on up- and down-regulated genes ----
# Universe: the tested genes. With 0 DEGs, run_go_enrichment() skips the
# analysis with a message.
genes_up   <- rownames(res_D_symbols)[res_D_symbols$diffexpressed == "up"]
genes_down <- rownames(res_D_symbols)[res_D_symbols$diffexpressed == "down"]
length(genes_up); length(genes_down)
gene_list_fc <- setNames(res_D_symbols$log2FoldChange, rownames(res_D_symbols))
universe_genes <- rownames(res_D_symbols)[!is.na(res_D_symbols$padj)]

go_down <- run_go_enrichment(genes_down, "down", plot_dir, gene_list_fc, label = cfg$group_label,
                             axis_prefix = cfg$axis_prefix, universe = universe_genes)
go_up   <- run_go_enrichment(genes_up, "up", plot_dir, gene_list_fc, label = cfg$group_label,
                             axis_prefix = cfg$axis_prefix, universe = universe_genes)
if (!is.null(go_up$go_bp))   head(go_up$go_bp@result)   else message("UP: no DEGs, enrichGO not run")
if (!is.null(go_down$go_bp)) head(go_down$go_bp@result) else message("DOWN: no DEGs, enrichGO not run")

# Functional GSEA (preranked) on all tested genes ----
# Hallmark and GO:BP, same function for all axes (run_functional_gsea()).
# Genes are ranked by the limma moderated t (column "stat").
gsea_fun <- run_functional_gsea(res_D_symbols, axis_prefix = cfg$axis_prefix,
                                tables_dir = tables_dir, plot_dir = plot_dir,
                                data_dir = data_dir, stat_col = "stat")
# Saves tables/D02_GSEA_*.csv, plots/D02_GSEA_*.png and data/D02_GSEA_functional.rds
gsea_fun$Hallmark[, c("pathway", "NES", "pval", "padj", "size")]
head(gsea_fun$GOBP[, c("pathway", "NES", "pval", "padj", "size")], 15)

# Objects kept in memory: go_up, go_down, gsea_fun.
