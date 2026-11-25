################################################################################
#                 SPECIFICITY OF THE AXIS A SIGNATURE TRANSFER
################################################################################
#namefile: Transferability/05_specificity.R
#creation dt: 14/10/2026
#
# Is the Axis A signature found in C (and D) because of the exposure, or would
# it be found anyway, because of how the signature is made or because of the
# correlation structure of blood data?
#
# The fgsea p-value is already a test against RANDOM gene sets of the same size
# (gene permutation). This test has two known weak points, checked here:
#
#   PART 1 - Random signatures matched for expression level.
#     The Axis A DEGs are not random genes: they tend to be more expressed than
#     average, and genes highly expressed in blood (ribosomal, mitochondrial,
#     histones, hemoglobins) often move together in any comparison. 1,000
#     random signatures with the same size AND the same distribution of mean
#     expression in A (baseMean, by decile) are built; the enrichment score (ES)
#     of each in the target ranking is compared with the observed ES -> empirical p.
#
#   PART 2 - Permutation of the sample labels (phenotype permutation,
#     Subramanian et al. 2005).
#     Gene permutation ignores the correlation between genes: if the signature
#     genes are correlated, a set can look enriched by chance. Permuting the
#     sample labels instead (poor / non-poor in C; adversity / control in D)
#     keeps the gene correlation and gives the ES distribution expected if the
#     exposure played no role. 1,000 permutations; in C the labels are
#     permuted within each sex (sex stays in the model).
#     For speed, the ranking of each permutation is computed with limma
#     (moderated t) on the VST matrix (C) or the normalized matrix (D); for C
#     the observed ES with the same method is also reported, so the comparison
#     is consistent. ES values are compared directly because each signature has
#     the same size in all permutations; NES = ES / mean |ES| of the permuted
#     values with the same sign (definition of Subramanian et al. 2005).
#
#   PART 3 - Per-sample signature score (mean z of the signature genes),
#     compared between groups with a linear model: it measures the SIZE of the
#     shift at the sample level (robust to co-expression).
#
# An external blood dataset without adversity is not used as a negative
# control: the choice of the contrast would be arbitrary, while the label
# permutation gives exactly the "no effect of the exposure" case on the same
# samples.
#
# Run from the project root, AFTER 01_signature_transfer_gsea.R:
#   source("Transferability/05_specificity.R", echo = TRUE)
# Time: about 15-30 minutes (preprocessing of C and D + 1,000 limma permutations on C).
# Outputs: Transferability/results/specificity_part1/2/3_*.csv, figures/specificity_null.png/.pdf
################################################################################

source("00_config.R")
source("00_functions.R")
load_pipeline_libraries()
library(fgsea)
set.seed(18)
N_RAND <- 1000
N_PERM <- 1000
res_dir <- "Transferability/results"
fig_dir <- "Transferability/figures"

# classic GSEA enrichment score (weight 1, as in fgsea)
es_of <- function(stats, genes) {
  sel <- match(intersect(genes, names(stats)), names(stats))
  fgsea::calcGseaStat(stats, selectedStats = sel, gseaParam = 1)
}
emp_p <- function(obs, null) {             # one-sided, in the direction of the observed value
  if (obs < 0) (1 + sum(null <= obs)) / (1 + length(null)) else (1 + sum(null >= obs)) / (1 + length(null))
}
nes_vs_null <- function(obs, null) {
  if (obs < 0) obs / mean(abs(null[null < 0])) else obs / mean(null[null > 0])
}

# Axis A signature and observed rankings ----
dict  <- build_gene_dictionary("Transferability/Dictionary")
A <- readRDS("Transferability/data/A_DEG_full_results_entrez.rds")
A <- A[!is.na(A$entrez) & !duplicated(A$entrez), ]
res_A <- readRDS(file.path(DIMENSIONS$mass_trauma$output_dir, "data", "DEG_full_results.rds"))
sig_up   <- unique(na.omit(to_entrez(rownames(res_A)[res_A$diffexpressed == "up"],   dict)$entrez))
sig_down <- unique(na.omit(to_entrez(rownames(res_A)[res_A$diffexpressed == "down"], dict)$entrez))
signatures <- list(WTC_down = sig_down, WTC_up = sig_up, WTC_all = unique(c(sig_up, sig_down)))

ranked <- list(C = build_ranked_stats(readRDS("Transferability/data/C_DEG_full_results_entrez.rds")),
               D = build_ranked_stats(readRDS("Transferability/data/D_DEG_full_results_entrez.rds")))
obs <- do.call(rbind, lapply(names(ranked), function(t) do.call(rbind, lapply(names(signatures), function(s)
  data.frame(target = t, signature = s, size = sum(signatures[[s]] %in% names(ranked[[t]])),
             ES_obs = es_of(ranked[[t]], signatures[[s]]))))))
print(obs, digits = 3, row.names = FALSE)

# Part 1: random signatures matched for expression in A (baseMean deciles) ----
part1 <- list(); null1 <- list()
for (t in names(ranked)) {
  pool <- A[A$entrez %in% names(ranked[[t]]), c("entrez", "baseMean")]
  pool$decile <- cut(rank(pool$baseMean, ties.method = "first"), 10, labels = FALSE)
  for (s in names(signatures)) {
    sig <- intersect(signatures[[s]], pool$entrez)
    target_dec <- table(factor(pool$decile[match(sig, pool$entrez)], levels = 1:10))
    cand <- split(pool$entrez[!pool$entrez %in% unlist(signatures)], pool$decile[!pool$entrez %in% unlist(signatures)])
    es_null <- replicate(N_RAND, {
      g <- unlist(lapply(names(target_dec), function(d) sample(cand[[d]], target_dec[[d]])))
      es_of(ranked[[t]], g)
    })
    o <- obs$ES_obs[obs$target == t & obs$signature == s]
    null1[[paste(t, s)]] <- es_null
    part1[[paste(t, s)]] <- data.frame(target = t, signature = s, test = "expression-matched random signatures",
                                       ES_obs = o, null_frac_same_sign = mean(sign(es_null) == sign(o)),
                                       null_q025 = quantile(es_null, .025), null_q975 = quantile(es_null, .975),
                                       NES_vs_null = nes_vs_null(o, es_null), p_empirical = emp_p(o, es_null))
  }
}
part1 <- do.call(rbind, part1); rownames(part1) <- NULL
print(part1, digits = 3, row.names = FALSE)
write.csv(part1, file.path(res_dir, "specificity_part1_matched_random.csv"), row.names = FALSE)

# Part 2: permutation of the sample labels (C and D) ----
run_preproc <- function(script) {
  env <- new.env(parent = globalenv()); sys.source(script, envir = env, chdir = FALSE); env
}
max_abs_by <- function(stat, key) {        # one statistic per key: the one with the largest |stat|
  ok <- !is.na(key) & !is.na(stat); stat <- stat[ok]; key <- key[ok]
  o <- order(-abs(stat)); stat <- stat[o]; key <- key[o]
  keep <- !duplicated(key); sort(setNames(stat[keep], key[keep]), decreasing = TRUE)
}
limma_t <- function(E, grp, covar = NULL) {
  d <- if (is.null(covar)) model.matrix(~ grp) else model.matrix(~ covar + grp)
  fit <- eBayes(lmFit(E, d)); fit$t[, ncol(d)]
}

## Axis C: VST, genes of the original IQR filter, limma t ~ sex + poverty ----
envC <- run_preproc("C_poverty_ses/00_preprocessing.R")
C_orig <- readRDS("Transferability/data/C_DEG_full_results_entrez.rds")
E_C <- envC$norm$vst[intersect(C_orig$original_id, rownames(envC$norm$vst)), ]
key_C <- setNames(C_orig$entrez, C_orig$original_id)[rownames(E_C)]
mC <- envC$norm$meta[colnames(E_C), ]
grp_C <- factor(mC[[DIMENSIONS$poverty_ses$condition_col]], levels = c("No", "Yes"))
sex_C <- factor(mC$sex)

## Axis D: normalized matrix, probes of the IQR filter, probe -> symbol -> Entrez ----
envD <- run_preproc("D_early_adversity/00_preprocessing.R")
cfgD <- DIMENSIONS$early_adversity
E_D_all <- envD$norm$exprs
iqr_D <- apply(E_D_all, 1, IQR); E_D <- E_D_all[iqr_D > quantile(iqr_D, cfgD$prc_iqr), ]
gplD <- GEOquery::Table(envD$gpl_D)
symc <- intersect(c("Gene symbol", "GENE_SYMBOL", "Symbol"), colnames(gplD))[1]
p2s <- setNames(trimws(sub("///.*$", "", as.character(gplD[[symc]]))), as.character(gplD$ID))
sym_D <- p2s[rownames(E_D)]; sym_D[sym_D == ""] <- NA
sym2ent <- to_entrez(unique(na.omit(sym_D)), dict)
sym2ent <- setNames(sym2ent$entrez, unique(na.omit(sym_D)))
mD <- envD$norm$meta[colnames(E_D), ]
grp_D <- factor(mD[[cfgD$condition_col]], levels = c("No", "Yes"))

rank_C <- function(g) max_abs_by(limma_t(E_C, g, sex_C), key_C)
rank_D <- function(g) {                     # probe -> symbol (max |t|) -> Entrez (max |t|)
  t_probe <- limma_t(E_D, g); s <- max_abs_by(t_probe, sym_D)
  max_abs_by(s, sym2ent[names(s)])
}
perm_within <- function(g, strata = NULL) {
  if (is.null(strata)) return(sample(g))
  out <- g; for (l in levels(strata)) { i <- which(strata == l); out[i] <- g[i][sample(length(i))] }; out
}

obs_limma <- list(C = rank_C(grp_C), D = rank_D(grp_D))
cat("Check D: correlation between the rebuilt and the original ranking =",
    round(cor(obs_limma$D[names(ranked$D)], ranked$D, use = "complete.obs"), 4), "(expected ~1)\n")
cat("C: Spearman correlation, limma vs DESeq2 =",
    round(cor(obs_limma$C[names(ranked$C)], ranked$C, use = "complete.obs", method = "spearman"), 3), "\n")

perm_es <- list()
for (t in c("C", "D")) {
  g0 <- if (t == "C") grp_C else grp_D; st <- if (t == "C") sex_C else NULL
  f  <- if (t == "C") rank_C else rank_D
  cat("Permutations, axis", t, "...\n")
  perm_es[[t]] <- t(replicate(N_PERM, {
    r <- f(perm_within(g0, st)); sapply(signatures, function(s) es_of(r, s))
  }))
}
part2 <- do.call(rbind, lapply(c("C", "D"), function(t) do.call(rbind, lapply(names(signatures), function(s) {
  o <- es_of(obs_limma[[t]], signatures[[s]]); nul <- perm_es[[t]][, s]
  data.frame(target = t, signature = s, test = "sample-label permutation",
             ES_obs = o, null_frac_same_sign = mean(sign(nul) == sign(o)),
             null_q025 = quantile(nul, .025), null_q975 = quantile(nul, .975),
             NES_vs_null = nes_vs_null(o, nul), p_empirical = emp_p(o, nul))
}))))
rownames(part2) <- NULL
print(part2, digits = 3, row.names = FALSE)
write.csv(part2, file.path(res_dir, "specificity_part2_label_permutation.csv"), row.names = FALSE)

# Part 3: per-sample signature score ----
# The ES measures WHERE the signature genes end up in the ranking, not HOW MUCH
# they shift. Here each person gets a score = mean standardized expression (z
# per gene) of the signature genes, and the difference between exposed and
# unexposed is tested with a linear model (C: adjusted for sex). The test is at
# the sample level, so it is robust to co-expression between genes. The power
# of a two-group t-test for the observed effect is also reported.
score_test <- function(E, key, grp, covar = NULL, label) {
  z <- t(scale(t(E))); ent <- key[rownames(z)]
  do.call(rbind, lapply(names(signatures), function(s) {
    sc <- colMeans(z[ent %in% signatures[[s]], , drop = FALSE])
    f <- if (is.null(covar)) lm(sc ~ grp) else lm(sc ~ covar + grp)
    cf <- summary(f)$coefficients[paste0("grp", levels(grp)[2]), ]
    d <- cf[1] / sd(sc); n1 <- sum(grp == levels(grp)[1]); n2 <- sum(grp == levels(grp)[2])
    pw <- pnorm(abs(d) / sqrt(1 / n1 + 1 / n2) - qnorm(.975))   # power (normal approximation)
    data.frame(target = label, signature = s, n_features = sum(ent %in% signatures[[s]]),   # D: probes (several probes per gene)
               diff_exposed_minus_unexposed = cf[1], score_sd = sd(sc), effect_size_SD = d,
               p = cf[4], power_two_sided_0.05 = pw)
  }))
}
sym2ent_probe <- sym2ent[sym_D]; names(sym2ent_probe) <- rownames(E_D)
part3 <- rbind(score_test(E_C, key_C, grp_C, sex_C, "C"),
               score_test(E_D, sym2ent_probe, grp_D, NULL, "D"))
rownames(part3) <- NULL
print(part3, digits = 3, row.names = FALSE)
write.csv(part3, file.path(res_dir, "specificity_part3_signature_score.csv"), row.names = FALSE)

# Figure: null distributions and observed ES (C and D, down and up components) ----
nl <- rbind(
  do.call(rbind, lapply(names(null1), function(k) { p <- strsplit(k, " ")[[1]]
    data.frame(target = p[1], signature = p[2], test = "Expression-matched random signatures", ES = null1[[k]]) })),
  do.call(rbind, lapply(c("C", "D"), function(t) do.call(rbind, lapply(names(signatures), function(s)
    data.frame(target = t, signature = s, test = "Sample-label permutation", ES = perm_es[[t]][, s]))))))
ob <- rbind(data.frame(part1[, c("target", "signature")], test = "Expression-matched random signatures", ES = part1$ES_obs),
            data.frame(part2[, c("target", "signature")], test = "Sample-label permutation", ES = part2$ES_obs))
keep <- function(d) d[d$signature != "WTC_all", ]
nl <- keep(nl); ob <- keep(ob)
lab <- c(C = "Axis C (poverty)", D = "Axis D (early adversity)")
nl$panel <- paste(lab[nl$target], sub("WTC_", "", nl$signature), sep = " - ")
ob$panel <- paste(lab[ob$target], sub("WTC_", "", ob$signature), sep = " - ")
g <- ggplot(nl, aes(ES)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white", linewidth = 0.2) +
  geom_vline(data = ob, aes(xintercept = ES), colour = "#2a78d6", linewidth = 0.8) +
  facet_grid(test ~ panel, scales = "free") +
  labs(x = "Enrichment score (ES)", y = "Count",
       title = "Specificity of the Axis A signature transfer",
       caption = "Grey: null distribution (1,000 draws). Blue line: observed ES.") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "specificity_null.png"), g, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "specificity_null.pdf"), g, width = 10, height = 5.5)
