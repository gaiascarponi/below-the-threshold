################################################################################
#         FUNCTIONAL ENRICHMENT (GO / enrichR) - MASS TRAUMA (Axis A)
################################################################################
#namefile: 02_enrichment.R
#creation dt: 18/08/2026
#
# Axis A, step 3 of 3: over-representation analysis (GO, Enrichr) of the
# up- and down-regulated genes and functional GSEA on all tested genes. Run
# right after 00_preprocessing.R and 01_DEG.R in the same R session (it
# reuses deg, cfg and the output folders).

# ORA (GO / enrichR) on up- and down-regulated genes ----
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
# Plots saved with the 02_ prefix
head(go_up$go_bp@result)
head(go_down$go_bp@result)

## GO Molecular Function, same design as BP and CC ----
library(clusterProfiler); library(org.Hs.eg.db)
run_mf <- function(genes, tag) {
  mf <- enrichGO(gene = unique(genes), OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 universe = union(unique(universe_genes), unique(genes)))
  write.csv(as.data.frame(mf@result),
            file.path("A_mass_trauma/03_results/tables", paste0("A02_GO_", tag, "_MF_enrichGO.csv")),
            row.names = FALSE)
  if (nrow(as.data.frame(mf)) > 0) {
    png(file.path(plot_dir, paste0("A02_GO_", tag, "_MF.png")), width = 1000, height = 800, res = 120)
    print(dotplot(mf, showCategory = 10) + ggtitle(paste("Molecular Function:", tag, "-regulated")))
    dev.off()
  }
  cat(tag, ": ", nrow(as.data.frame(mf)), " MF terms with padj < 0.05\n")
  mf
}
mf_down <- run_mf(genes_down, "DOWN")
mf_up   <- run_mf(genes_up,   "UP")

## Check of GO:CC of the down-regulated genes, whole-genome vs tested-genes background ----
library(clusterProfiler); library(org.Hs.eg.db)
cc_genome <- enrichGO(gene = unique(genes_down), OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "CC", pAdjustMethod = "BH")           # no universe
cc_tested <- enrichGO(gene = unique(genes_down), OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "CC", pAdjustMethod = "BH",
                      universe = union(unique(universe_genes), unique(genes_down)))
cmp <- function(x, lab) {
  r <- as.data.frame(x@result); r$rank <- seq_len(nrow(r))
  r <- r[grepl("ribosom", r$Description), c("Description", "GeneRatio", "BgRatio", "p.adjust", "rank")]
  cat("\n==", lab, "==\n"); print(head(r, 6), row.names = FALSE)
}
cmp(cc_genome, "whole genome"); cmp(cc_tested, "tested genes")
head(as.data.frame(cc_genome)[, c("Description", "GeneRatio", "p.adjust")], 10)   # top 10 with the whole-genome background

## Check of GO:BP of the down-regulated genes, whole-genome vs tested-genes background ----
library(clusterProfiler); library(org.Hs.eg.db)
cc_genome <- enrichGO(gene = unique(genes_down), OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "BP", pAdjustMethod = "BH")           # no universe
cc_tested <- enrichGO(gene = unique(genes_down), OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "BP", pAdjustMethod = "BH",
                      universe = union(unique(universe_genes), unique(genes_down)))
cmp <- function(x, lab) {
  r <- as.data.frame(x@result); r$rank <- seq_len(nrow(r))
  r <- r[grepl("ribosom", r$Description), c("Description", "GeneRatio", "BgRatio", "p.adjust", "rank")]
  cat("\n==", lab, "==\n"); print(head(r, 6), row.names = FALSE)
}
cmp(cc_genome, "whole genome"); cmp(cc_tested, "tested genes")
head(as.data.frame(cc_genome)[, c("Description", "GeneRatio", "p.adjust")], 10)   # top 10 with the whole-genome background

# Functional GSEA (preranked) on all tested genes ----
# Hallmark and GO:BP, same function for all axes (run_functional_gsea()).
gsea_fun <- run_functional_gsea(deg$res_filtered, axis_prefix = cfg$axis_prefix,
                                tables_dir = tables_dir, plot_dir = plot_dir,
                                data_dir = data_dir, stat_col = "stat")
gsea_fun$Hallmark[, c("pathway", "NES", "pval", "padj", "size")]
head(gsea_fun$GOBP[, c("pathway", "NES", "pval", "padj", "size")], 15)
# Saves tables/<P>02_GSEA_*.csv, plots/<P>02_GSEA_*.png, data/<P>02_GSEA_functional.rds.

# Objects kept in memory: go_up, go_down, gsea_fun (plus raw, metadata,
# pre, norm and deg from the previous scripts).
