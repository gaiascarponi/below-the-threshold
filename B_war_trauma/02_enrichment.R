################################################################################
#      FUNCTIONAL ENRICHMENT (GO / enrichR) - WAR/COMBAT TRAUMA (Axis B)
################################################################################
#namefile: 02_enrichment.R
#creation dt: 18/08/2026
#
# Axis B, step 3 of 3: over-representation analysis and functional GSEA.
# Run right after 01_DEG.R in the same R session (it reuses deg, cfg and
# the output folders).

# ORA (GO / enrichR) on up- and down-regulated genes ----
# Gene identifier type of this dataset (symbols or Ensembl IDs)
head(rownames(deg$res_filtered))

genes_up   <- deg$sig_genes[deg$res_filtered[deg$sig_genes, "diffexpressed"] == "up"]
genes_down <- deg$sig_genes[deg$res_filtered[deg$sig_genes, "diffexpressed"] == "down"]
length(genes_up)
length(genes_down)

go_up <- run_go_enrichment(genes_up, direction = "up", plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)
go_down <- run_go_enrichment(genes_down, direction = "down", plot_dir = plot_dir, axis_prefix = cfg$axis_prefix)

# With 0 DEGs (Axis B) the ORA cannot run: run_go_enrichment() returns NULL
# and a message is shown. The functional analysis of this axis is the
# preranked GSEA below, which needs no significance threshold.
if (!is.null(go_up$go_bp))   go_up$go_bp@result   else message("UP: no DEGs, enrichGO not run")
if (!is.null(go_down$go_bp)) go_down$go_bp@result else message("DOWN: no DEGs, enrichGO not run")

# Functional GSEA (preranked) on all tested genes ----
# Hallmark and GO:BP, same function for all axes (run_functional_gsea()).
gsea_fun <- run_functional_gsea(deg$res_filtered, axis_prefix = cfg$axis_prefix,
                                tables_dir = tables_dir, plot_dir = plot_dir,
                                data_dir = data_dir, stat_col = "stat")
# Saves tables/<P>02_GSEA_*.csv, plots/<P>02_GSEA_*.png and data/<P>02_GSEA_functional.rds
gsea_fun$Hallmark[, c("pathway", "NES", "pval", "padj", "size")]
head(gsea_fun$GOBP[, c("pathway", "NES", "pval", "padj", "size")], 15)

# Objects kept in memory: go_up, go_down, gsea_fun.
