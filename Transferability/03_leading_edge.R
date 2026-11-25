################################################################################
#          LEADING EDGE OF THE AXIS A SIGNATURE TRANSFER (A->C, A->D)
################################################################################
#namefile: Transferability/03_leading_edge.R
#creation dt: 14/10/2026
#
# Which genes does Axis A share with C (and D)? The leading edge of a GSEA is
# the subset of signature genes found before the peak of the running sum: the
# genes that drive the enrichment (Subramanian et al. 2005).
#
# For each significant leading edge (A->C down, A->C up, A->D down):
#   1. gene table: symbol, log2FC and padj in A, statistic, log2FC, nominal p
#      and rank in C (and D, B), switch gene in A, degree in the STRING network
#      of A, function (STRING / UniProt annotation)
#   2. shared core: genes in the leading edge of both A->C and A->D (down)
#   3. functional enrichment of the leading edge (ORA, hypergeometric test) on
#      MSigDB Hallmark and GO:BP; universe = Axis A signature genes present in
#      the target ranking (which functions separate the genes that transfer
#      from the signature genes that do not)
#   4. switch genes: are they over-represented in the leading edge? (Fisher)
#   5. figure: log2FC in A vs statistic in C for the signature genes, with the
#      leading edge highlighted
#
# Requirements: 00_harmonize_gene_ids.R and 01_signature_transfer_gsea.R
# already run; Axis A network results in A_mass_trauma/04_network/ (switch
# genes, STRING).
#
# Run from the project root:
#   source("Transferability/03_leading_edge.R", echo = TRUE)
# Outputs: Transferability/results/leading_edge_*.csv,
#          Transferability/figures/leading_edge_A_to_C.png/.pdf
################################################################################

library(ggplot2)
set.seed(18)

res_dir <- "Transferability/results"
fig_dir <- "Transferability/figures"
net_dir <- "A_mass_trauma/04_network"

# Inputs ----
A <- readRDS("Transferability/data/A_DEG_full_results_entrez.rds")
C <- readRDS("Transferability/data/C_DEG_full_results_entrez.rds")
D <- readRDS("Transferability/data/D_DEG_full_results_entrez.rds")
B <- readRDS("Transferability/data/B_DEG_full_results_entrez.rds")
gsea <- readRDS(file.path(res_dir, "gsea_signature_transfer.rds"))

switch_genes <- readLines(file.path(net_dir, "input", "A_switch_genes.txt"))
switch_genes <- trimws(switch_genes[nzchar(trimws(switch_genes))])
deg_str <- read.delim(file.path(net_dir, "string", "DEG_all", "string_node_degrees_418_nodes.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
ann_str <- read.delim(file.path(net_dir, "string", "DEG_all", "string_protein_annotations_418_nodes.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE, quote = "")
names(deg_str)[1] <- "node"; names(ann_str)[1] <- "node"
cat("Switch genes:", length(switch_genes), "| nodi STRING:", nrow(deg_str), "\n")

# ranking (statistic) and rank in the target; 1 = top
rank_tab <- function(R, prefix) {
  R <- R[!is.na(R$entrez) & !is.na(R$stat) & !duplicated(R$entrez), ]
  R$rank <- rank(-R$stat, ties.method = "first")
  pv <- if ("pvalue" %in% names(R)) R$pvalue else rep(NA_real_, nrow(R))   # D (limma): padj only
  out <- data.frame(entrez = R$entrez, stat = R$stat, log2FC = R$log2FoldChange,
                    pvalue = pv, rank = R$rank, n = nrow(R))
  names(out)[-1] <- paste0(prefix, "_", names(out)[-1])
  out
}
rk <- Reduce(function(x, y) merge(x, y, by = "entrez", all.x = TRUE),
             list(rank_tab(C, "C"), rank_tab(D, "D"), rank_tab(B, "B")))
rk <- merge(rank_tab(C, "C")[, "entrez", drop = FALSE], rk, by = "entrez", all = TRUE)

A_sig <- A[A$diffexpressed %in% c("up", "down") & !is.na(A$entrez), ]
A_sig <- A_sig[!duplicated(A_sig$entrez), ]
cat("Axis A signature (Entrez):", table(A_sig$diffexpressed), "\n")

# short annotation: first sentence of the STRING description (protein name).
# In the STRING file the headers of columns 3 and 4 are swapped with respect
# to the data: the text column without URLs is used.
ann_col <- names(ann_str)[3:4][which.min(sapply(ann_str[3:4], function(x) mean(grepl("^https?://", x))))]
short_ann <- function(sym) {
  a <- ann_str[[ann_col]][match(sym, ann_str$node)]
  a <- sub(";.*$", "", a)
  ifelse(is.na(a), "", a)
}

# Leading-edge tables ----
le_genes <- function(target, pathway) {
  g <- gsea[[target]]; strsplit(g$leadingEdge[g$pathway == pathway], ";")[[1]]
}
make_le_table <- function(target, pathway) {
  le <- le_genes(target, pathway)
  a <- A_sig[match(le, A_sig$entrez), ]
  t <- data.frame(entrez = le, symbol = a$symbol_current,
                  A_direction = a$diffexpressed, A_log2FC = a$log2FoldChange, A_padj = a$padj)
  t <- merge(t, rk, by = "entrez", all.x = TRUE, sort = FALSE)
  t$switch_gene  <- t$symbol %in% switch_genes
  t$string_degree <- deg_str$node_degree[match(t$symbol, deg_str$node)]
  t$function_string <- short_ann(t$symbol)
  t$histone_hb <- grepl("^H[1-4]|^H2A|^H2B|^H3|^H4|^HB[ABDEGMQZ]", t$symbol)
  ord <- if (grepl("down", pathway)) order(t[[paste0(sub("A_to_", "", target), "_stat")]])
         else order(-t[[paste0(sub("A_to_", "", target), "_stat")]])
  t[ord, ]
}
le_C_down <- make_le_table("A_to_C", "WTC_down")
le_C_up   <- make_le_table("A_to_C", "WTC_up")
le_D_down <- make_le_table("A_to_D", "WTC_down")
for (nm in c("le_C_down", "le_C_up", "le_D_down")) {
  t <- get(nm)
  cat("\n", nm, ": ", nrow(t), " genes | switch genes: ", sum(t$switch_gene),
      " | histones/hemoglobins: ", sum(t$histone_hb), "\n", sep = "")
  print(head(t[, c("symbol", "A_log2FC", "C_stat", "C_rank", "D_stat", "switch_gene",
                   "string_degree", "function_string")], 15), digits = 3, row.names = FALSE)
}
write.csv(le_C_down, file.path(res_dir, "leading_edge_A_to_C_down.csv"), row.names = FALSE)
write.csv(le_C_up,   file.path(res_dir, "leading_edge_A_to_C_up.csv"),   row.names = FALSE)
write.csv(le_D_down, file.path(res_dir, "leading_edge_A_to_D_down.csv"), row.names = FALSE)

# Shared core (leading edges A->C down and A->D down) ----
core <- le_C_down[le_C_down$entrez %in% le_D_down$entrez, ]
cat("\nShared core A->C and A->D (down):", nrow(core), "genes out of", nrow(le_C_down), "e", nrow(le_D_down), "\n")
print(core[, c("symbol", "A_log2FC", "C_stat", "D_stat", "switch_gene", "function_string")],
      digits = 3, row.names = FALSE)
# expected by chance: hypergeometric on the down-regulated signature genes present in both rankings
down_both <- A_sig$entrez[A_sig$diffexpressed == "down" & A_sig$entrez %in% C$entrez & A_sig$entrez %in% D$entrez]
k <- nrow(core); n1 <- sum(le_C_down$entrez %in% down_both); n2 <- sum(le_D_down$entrez %in% down_both)
p_core <- phyper(k - 1, n1, length(down_both) - n1, n2, lower.tail = FALSE)
cat("Expected by chance:", round(n1 * n2 / length(down_both), 1), "| observed:", k,
    "| hypergeometric p =", signif(p_core, 3), "\n")
write.csv(core, file.path(res_dir, "leading_edge_core_C_D_down.csv"), row.names = FALSE)

# Functional enrichment of the leading edge (ORA) ----
H  <- readRDS("Transferability/Dictionary/msigdb_H.rds")
BP <- readRDS("Transferability/Dictionary/msigdb_C5_GO_BP.rds")
ora <- function(genes, universe, sets, min_size = 5) {
  sets <- lapply(sets, intersect, universe)
  sets <- sets[lengths(sets) >= min_size]
  r <- do.call(rbind, lapply(names(sets), function(s) {
    k <- length(intersect(genes, sets[[s]])); m <- length(sets[[s]])
    data.frame(set = s, k = k, set_size = m, expected = length(genes) * m / length(universe),
               p = phyper(k - 1, m, length(universe) - m, length(genes), lower.tail = FALSE))
  }))
  r$padj <- p.adjust(r$p, "BH"); r[order(r$p), ]
}
ora_runs <- list(
  C_down = list(le = le_C_down$entrez, uni = intersect(A_sig$entrez[A_sig$diffexpressed == "down"], C$entrez)),
  C_up   = list(le = le_C_up$entrez,   uni = intersect(A_sig$entrez[A_sig$diffexpressed == "up"],   C$entrez)),
  D_down = list(le = le_D_down$entrez, uni = intersect(A_sig$entrez[A_sig$diffexpressed == "down"], D$entrez)))
# two universes: (i) Axis A signature genes present in the target -> what
# separates the genes that transfer from the rest of the signature; (ii) all
# genes tested in the target -> what the leading edge is about, compared with
# the expressed genome
uni_all <- list(C_down = C$entrez, C_up = C$entrez, D_down = D$entrez)
ora_out <- do.call(rbind, lapply(names(ora_runs), function(nm) {
  x <- ora_runs[[nm]]; ua <- unique(na.omit(uni_all[[nm]]))
  rbind(data.frame(leading_edge = nm, universe = "signature", collection = "Hallmark", ora(x$le, x$uni, H)),
        data.frame(leading_edge = nm, universe = "signature", collection = "GO:BP",    ora(x$le, x$uni, BP)),
        data.frame(leading_edge = nm, universe = "all_tested", collection = "Hallmark", ora(x$le, ua, H)),
        data.frame(leading_edge = nm, universe = "all_tested", collection = "GO:BP",    ora(x$le, ua, BP)))
}))
cat("\nORA of the leading edge (universe = Axis A signature present in the target), top by p:\n")
for (nm in names(ora_runs)) for (u in c("signature", "all_tested")) {
  x <- ora_out[ora_out$leading_edge == nm & ora_out$universe == u, ]
  x <- x[order(x$p), ]
  cat("--", nm, "| universe:", u, "\n")
  print(head(x[, c("collection", "set", "k", "set_size", "expected", "p", "padj")], 8), digits = 3, row.names = FALSE)
}
write.csv(ora_out, file.path(res_dir, "leading_edge_ORA.csv"), row.names = FALSE)

# Switch genes in the leading edge (Fisher) ----
sw_test <- function(le_tab, direction, target_df) {
  uni <- A_sig[A_sig$diffexpressed == direction & A_sig$entrez %in% target_df$entrez, ]
  in_le <- uni$entrez %in% le_tab$entrez; is_sw <- uni$symbol_current %in% switch_genes
  ft <- fisher.test(table(factor(in_le, c(TRUE, FALSE)), factor(is_sw, c(TRUE, FALSE))))
  data.frame(switch_in_signature = sum(is_sw), switch_in_le = sum(in_le & is_sw),
             le_size = sum(in_le), signature_size = nrow(uni),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}
sw <- rbind(cbind(leading_edge = "C_down", sw_test(le_C_down, "down", C)),
            cbind(leading_edge = "C_up",   sw_test(le_C_up,   "up",   C)),
            cbind(leading_edge = "D_down", sw_test(le_D_down, "down", D)))
print(sw, digits = 3, row.names = FALSE)
write.csv(sw, file.path(res_dir, "leading_edge_switch_genes.csv"), row.names = FALSE)

# Figure: agreement between A and C for the signature genes ----
pd <- merge(A_sig[, c("entrez", "symbol_current", "diffexpressed", "log2FoldChange")],
            rank_tab(C, "C")[, c("entrez", "C_stat")], by = "entrez")
pd$group <- ifelse(pd$entrez %in% le_C_down$entrez, "Leading edge (down)",
            ifelse(pd$entrez %in% le_C_up$entrez, "Leading edge (up)", "Rest of the signature"))
pd$group <- factor(pd$group, levels = c("Rest of the signature", "Leading edge (down)", "Leading edge (up)"))
pd <- pd[order(pd$group), ]
lab <- rbind(head(le_C_down, 6), head(le_C_up, 6))
lab <- merge(lab[, c("entrez", "symbol")], pd, by = "entrez")
COL <- c("Rest of the signature" = "grey75", "Leading edge (down)" = "#2a78d6", "Leading edge (up)" = "#eb6834")
SHP <- c("Rest of the signature" = 16, "Leading edge (down)" = 16, "Leading edge (up)" = 17)
g <- ggplot(pd, aes(log2FoldChange, C_stat, colour = group, shape = group)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.85) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = lab, aes(label = symbol), size = 2.8, colour = "grey15",
                               show.legend = FALSE, min.segment.length = 0, seed = 18)
    else geom_text(data = lab, aes(label = symbol), size = 2.8, vjust = -0.8, show.legend = FALSE, colour = "grey15") } +
  scale_colour_manual(values = COL, name = NULL) + scale_shape_manual(values = SHP, name = NULL) +
  labs(x = "log2 fold change in Axis A (WTC, PTSD vs controls)",
       y = "Wald statistic in Axis C (poverty vs no poverty)",
       title = "Genes of the Axis A signature in the Axis C ranking",
       caption = "Leading edge of the A->C transfer (fgsea). Labels: top 6 leading-edge genes per direction.") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "leading_edge_A_to_C.png"), g, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "leading_edge_A_to_C.pdf"), g, width = 7.5, height = 6)
