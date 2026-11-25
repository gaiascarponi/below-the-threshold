################################################################################
#        CROSS-AXIS HEATMAP OF THE HALLMARK SETS (FUNCTIONAL GSEA NES)
################################################################################
#namefile: plot_hallmark_heatmap.R
#creation dt: 19/10/2026 (R version of plot_hallmark_heatmap.py)
#
# Reads the Hallmark GSEA tables of the 4 axes (<axis>02_GSEA_Hallmark.csv),
# keeps the sets significant (padj < 0.05) in at least one axis, orders them by
# hierarchical clustering (average linkage, Euclidean distance on the NES) and
# draws the NES heatmap with asterisks for significance.
#
# Run from the project root (open Project.Rproj):
#   source("Transferability/plot_hallmark_heatmap.R")
# Outputs:
#   Transferability/results/hallmark_NES_matrix.csv
#   Transferability/figures/X_GSEA_Hallmark_heatmap.png / .pdf
#
# The row order can differ slightly from the Python version (scipy and
# hclust() order tied leaves differently); the matrix values are identical.
################################################################################
library(ggplot2)

base_dir <- "."
out_dir  <- file.path(base_dir, "Transferability")

axes <- data.frame(
  ax    = c("A", "B", "C", "D"),
  dir   = c("A_mass_trauma", "B_war_trauma", "C_poverty_ses", "D_early_adversity"),
  label = c("A\nMass trauma\n(WTC, PTSD)", "B\nWar trauma\n(Marines, PTSD)",
            "C\nPoverty\n(HANDLS)", "D\nEarly adversity\n(monocytes)"),
  stringsAsFactors = FALSE)

# Readable names of the Hallmark sets (same replacements as the Python version)
pretty_hallmark <- function(x) {
  k <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x)))
  k <- paste0(toupper(substr(k, 1, 1)), substring(k, 2))
  rep <- c("Tnfa signaling via nfkb" = "TNFα signaling via NF-κB",
           "Il6 jak stat3" = "IL6/JAK/STAT3", "Il2 stat5" = "IL2/STAT5",
           "Myc " = "MYC ", "E2f " = "E2F ", "G2m " = "G2M ", "Mtorc1" = "mTORC1",
           "Pi3k akt mtor" = "PI3K/AKT/mTOR", "Kras " = "KRAS ",
           "Tgf beta" = "TGF-β", "Wnt beta catenin" = "WNT/β-catenin",
           "Uv response" = "UV response", "Dna repair" = "DNA repair",
           "P53 " = "p53 ", " dn" = " down")
  for (i in seq_along(rep)) k <- gsub(names(rep)[i], rep[i], k, fixed = TRUE)
  k
}

# NES and padj of each axis ----
long <- do.call(rbind, lapply(seq_len(nrow(axes)), function(i) {
  f <- file.path(base_dir, axes$dir[i], "03_results", "tables",
                 paste0(axes$ax[i], "02_GSEA_Hallmark.csv"))
  r <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(hallmark = pretty_hallmark(r$pathway), ax = axes$ax[i],
             NES = as.numeric(r$NES),
             padj = ifelse(is.na(suppressWarnings(as.numeric(r$padj))), 1,
                           as.numeric(r$padj)),
             stringsAsFactors = FALSE)
}))

M <- tapply(long$NES,  list(long$hallmark, long$ax), function(v) v[1])
P <- tapply(long$padj, list(long$hallmark, long$ax), function(v) v[1])
M <- M[order(rownames(M)), axes$ax, drop = FALSE]
P <- P[rownames(M), axes$ax, drop = FALSE]
P[is.na(P)] <- 1

# Sets significant (padj < 0.05) in at least one axis ----
keep <- apply(P < 0.05, 1, any)
M <- M[keep, , drop = FALSE]; P <- P[keep, , drop = FALSE]

# Order by hierarchical clustering (average, Euclidean) ----
M0 <- M; M0[is.na(M0)] <- 0
ord <- hclust(dist(M0, method = "euclidean"), method = "average")$order
M <- M[ord, , drop = FALSE]; P <- P[ord, , drop = FALSE]
cat(nrow(M), "Hallmark sets plotted\n")

# Save the matrix ----
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
mat_out <- data.frame(hallmark = rownames(M),
                      setNames(as.data.frame(round(M, 3)), paste0("NES_", axes$ax)),
                      setNames(as.data.frame(signif(P, 3)), paste0("padj_", axes$ax)),
                      check.names = FALSE)
# written explicitly in UTF-8 (names contain Greek letters), so it also works
# when the R session has no UTF-8 locale
csv_lines <- capture.output(write.csv(mat_out, row.names = FALSE))
con <- file(file.path(out_dir, "results", "hallmark_NES_matrix.csv"), open = "w", encoding = "UTF-8")
writeLines(enc2utf8(csv_lines), con)
close(con)

# Heatmap ----
vmax <- max(abs(M), na.rm = TRUE)
df <- data.frame(hallmark = rep(rownames(M), times = ncol(M)),
                 ax       = rep(colnames(M), each = nrow(M)),
                 NES      = as.vector(M),
                 padj     = as.vector(P),
                 stringsAsFactors = FALSE)
df$hallmark <- factor(df$hallmark, levels = rev(rownames(M)))   # first row at the top
df$ax       <- factor(df$ax, levels = axes$ax, labels = axes$label)
df$star     <- ifelse(df$padj < 0.001, "**", ifelse(df$padj < 0.05, "*", ""))
df$star_col <- ifelse(!is.na(df$NES) & abs(df$NES) > 0.62 * vmax, "white", "#212121")

p <- ggplot(df, aes(x = ax, y = hallmark, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = star, colour = star_col), size = 3.2, vjust = 0.75) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#2E6DA4", mid = "#F2F2F2", high = "#C0392B",
                       midpoint = 0, limits = c(-vmax, vmax), na.value = "grey90",
                       name = "NES (red: higher\nin exposed group)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = "* adjusted p < 0.05   ** adjusted p < 0.001   (Hallmark sets significant in at least one axis)") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.caption = element_text(size = 7, colour = "#555555", hjust = 0),
        plot.caption.position = "plot")

h <- 0.27 * nrow(M) + 1.9
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "figures", "X_GSEA_Hallmark_heatmap.png"), p,
       width = 6.4, height = h, dpi = 300, bg = "white")
# PDF with cairo for the Greek letters (alpha, kappa, beta); if cairo is not
# available the standard pdf() device is used
pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
ggsave(file.path(out_dir, "figures", "X_GSEA_Hallmark_heatmap.pdf"), p,
       width = 6.4, height = h, device = pdf_dev)
