################################################################################
#                  SOCIAL COST MODEL (DALY): SHARED FUNCTIONS
################################################################################
#namefile: DALY/00_daly_functions.R
#creation dt: 13/10/2026
#
# Functions used by the axis scripts (01_axisA_..., 02_axisB_..., ...):
#   1. inputs (GBD disability weights, BLS CPI-U)
#   2. distributions for the Monte Carlo simulation (quantile matching)
#   3. DALY model of PTSD (Axes A and B: same disease, different populations)
#   4. output tables and plots
#   5. attributable fractions and GBD data by age (Axes C and D)
#
# Distributions (see Chapter 3, "Uncertainty"):
#   - disability weights and prevalences are proportions in [0, 1] -> BETA,
#     with (a, b) chosen by quantile matching: the 2.5%, 50% and 97.5%
#     quantiles must match the published estimate and interval
#   - severity proportions (4 categories summing to 1):
#       main        = independent betas (quantile matching) + renormalization
#       sensitivity = Dirichlet with a single concentration K
#     The Dirichlet has one dispersion parameter for all categories and
#     cannot reproduce intervals of very different width at the same time
#     (asymptomatic 27.5-30.1% vs mild 34.2-44.2%), so it is only a
#     sensitivity analysis.
#   - monetary value of a DALY: UNIFORM over the range
# For each input a table of published vs reproduced values (maximum
# relative error over the three quantiles) is saved.
#
# Use: source("DALY/00_daly_functions.R") at the top of each axis script
################################################################################

library(readxl)
library(ggplot2)

DALY_DIR <- "DALY"
SRC_DIR  <- file.path(DALY_DIR, "data", "source")
RES_DIR  <- file.path(DALY_DIR, "results")
FIG_DIR  <- file.path(DALY_DIR, "figures")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

DALY_VALUE <- c(low = 100000, high = 150000)   # USD per DALY (US range)
TARGET_YR  <- 2024                             # year of the output dollars

# Severity distribution of anxiety disorders, GBD 2019
# (Yang et al. 2021, Epidemiol Psychiatr Sci 30:e36, supplementary material)
SEVERITY_ANXIETY <- rbind(asymptomatic = c(est = .288, low = .275, high = .301),
                          mild         = c(est = .393, low = .342, high = .442),
                          moderate     = c(est = .191, low = .158, high = .227),
                          severe       = c(est = .127, low = .092, high = .167))

# Inputs ----

# GBD 2019 disability weights (IHME 2020): cells like "0.133\n(0.091-0.186)"
read_disability_weights <- function(sequelae,
                                    file = file.path(SRC_DIR, "IHME_GBD_2019_DISABILITY_WEIGHTS_Y2020M010D15.XLSX")) {
  raw <- suppressMessages(read_excel(file, skip = 1))
  names(raw)[1:4] <- c("sequela", "health_state", "lay_description", "dw")
  parse <- function(x) {
    n <- as.numeric(regmatches(x, gregexpr("[0-9]*\\.?[0-9]+", x))[[1]])
    setNames(n[1:3], c("est", "low", "high"))
  }
  out <- t(sapply(sequelae, function(s) {
    hit <- raw$dw[raw$sequela == s]
    if (length(hit) != 1) stop("Sequela not found or duplicated: ", s)
    parse(hit)
  }))
  rownames(out) <- names(sequelae)
  out
}

# CPI-U (BLS, series CUUR0000SA0): annual mean of the 12 monthly values;
# with several SeriesReport files the most recent one is used (timestamp in the name)
read_cpi <- function(file = tail(sort(list.files(SRC_DIR, "^SeriesReport.*\\.xlsx$", full.names = TRUE)), 1)) {
  raw <- suppressMessages(read_excel(file, col_names = FALSE))
  hdr <- which(raw[[1]] == "Year")
  x   <- raw[(hdr + 1):nrow(raw), 1:13]
  data.frame(year   = as.integer(as.numeric(x[[1]])),
             annual = rowMeans(sapply(x[, -1], function(v) suppressWarnings(as.numeric(v)))))
}

inflation_factor <- function(cpi, from_year, to_year = TARGET_YR) {
  f <- cpi$annual[cpi$year == to_year] / cpi$annual[cpi$year == from_year]
  if (length(f) != 1 || is.na(f)) stop("CPI missing for year ", from_year, " or ", to_year)
  f
}

# Distributions (quantile matching) ----

PROBS <- c(.025, .5, .975)

# maximum relative error between theoretical quantiles and published values (%)
rel_err <- function(q, target) 100 * max(abs(q - target) / target)

# BETA: (a, b) such that qbeta(2.5%, 50%, 97.5%) ~ (low, est, high)
fit_beta <- function(est, low, high) {
  est <- unname(est); low <- unname(low); high <- unname(high)   # drop names like "a.ace1"
  target <- c(low, est, high)
  obj <- function(p) {
    q <- qbeta(PROBS, exp(p[1]), exp(p[2]))
    sum(((q - target) / target)^2)
  }
  sd0 <- (high - low) / (2 * 1.96)                    # start: method of moments
  k0  <- max(est * (1 - est) / sd0^2 - 1, 2)
  fit <- optim(log(c(est * k0, (1 - est) * k0)), obj,
               control = list(reltol = 1e-14, maxit = 5000))
  ab <- exp(fit$par)
  c(a = ab[1], b = ab[2], err_pct = rel_err(qbeta(PROBS, ab[1], ab[2]), target))
}

# comparison with the truncated normal used in the previous version
truncnorm_err <- function(est, low, high) {
  sd <- (high - low) / (2 * 1.96)
  p0 <- pnorm(0, est, sd); p1 <- pnorm(1, est, sd)
  q  <- qnorm(p0 + PROBS * (p1 - p0), est, sd)
  rel_err(q, c(low, est, high))
}

# published vs reproduced check table for an est/low/high matrix
fit_table <- function(M, group) {
  do.call(rbind, lapply(rownames(M), function(r) {
    b <- fit_beta(M[r, "est"], M[r, "low"], M[r, "high"])
    q <- qbeta(PROBS, b["a"], b["b"])
    data.frame(group = group, parameter = r,
               published_est = M[r, "est"], published_low = M[r, "low"], published_high = M[r, "high"],
               beta_a = b["a"], beta_b = b["b"],
               beta_q025 = q[1], beta_q50 = q[2], beta_q975 = q[3],
               beta_err_pct = b["err_pct"],
               truncnorm_err_pct = truncnorm_err(M[r, "est"], M[r, "low"], M[r, "high"]),
               row.names = NULL)
  }))
}

# DIRICHLET: concentration K such that the marginals Beta(K*m, K*(1-m))
# reproduce the published intervals as closely as possible (sensitivity)
fit_dirichlet <- function(S) {
  m <- S[, "est"] / sum(S[, "est"])
  obj <- function(lk) {
    K <- exp(lk)
    sum(sapply(seq_along(m), function(i) {
      q <- qbeta(c(.025, .975), K * m[i], K * (1 - m[i]))
      sum(((q - S[i, c("low", "high")]) / S[i, c("low", "high")])^2)
    }))
  }
  K <- exp(optimize(obj, c(0, 12))$minimum)
  err <- sapply(seq_along(m), function(i)
    rel_err(qbeta(PROBS, K * m[i], K * (1 - m[i])), S[i, c("low", "est", "high")]))
  list(alpha = setNames(K * m, rownames(S)), K = K, err_pct = setNames(err, rownames(S)))
}

rdirichlet <- function(n, alpha) {
  g <- sapply(alpha, function(a) rgamma(n, shape = a))
  g / rowSums(g)
}

rbeta_fit <- function(n, fit) rbeta(n, fit["a"], fit["b"])

# DALY model of PTSD (Axes A and B) ----
#   YLD per case = sum_s p_s * DW_s   (asymptomatic weight 0; YLL = 0, as in GBD)
#   cases        = Pop_M * prev_M + Pop_F * prev_F
#   total DALYs  = YLD per case * cases
#   value ($)    = YLD per case * value of a DALY ~ U(100k, 150k)
run_ptsd_daly <- function(DW, SEV, PREV, pop_m, pop_f, n_sim = 10000, seed = 18,
                          severity_method = c("beta", "dirichlet")) {
  severity_method <- match.arg(severity_method)
  set.seed(seed)
  fdw  <- lapply(rownames(DW),   function(r) fit_beta(DW[r, "est"], DW[r, "low"], DW[r, "high"]));   names(fdw)  <- rownames(DW)
  fprv <- lapply(rownames(PREV), function(r) fit_beta(PREV[r, "est"], PREV[r, "low"], PREV[r, "high"])); names(fprv) <- rownames(PREV)

  dw_sim <- sapply(names(fdw), function(s) rbeta_fit(n_sim, fdw[[s]]))
  if (severity_method == "beta") {
    fsev    <- lapply(rownames(SEV), function(r) fit_beta(SEV[r, "est"], SEV[r, "low"], SEV[r, "high"]))
    sev_sim <- sapply(fsev, function(f) rbeta_fit(n_sim, f))
    colnames(sev_sim) <- rownames(SEV)
    sev_sim <- sev_sim / rowSums(sev_sim)
  } else {
    sev_sim <- rdirichlet(n_sim, fit_dirichlet(SEV)$alpha)
  }
  yld_case <- rowSums(sev_sim[, c("mild", "moderate", "severe")] * dw_sim[, c("mild", "moderate", "severe")])
  n_cases  <- pop_m * rbeta_fit(n_sim, fprv[["M"]]) + pop_f * rbeta_fit(n_sim, fprv[["F"]])
  data.frame(yld_case  = yld_case,
             n_cases   = n_cases,
             yld_total = yld_case * n_cases,
             rate_100k = yld_case * n_cases / (pop_m + pop_f) * 1e5,
             usd_case  = yld_case * runif(n_sim, DALY_VALUE["low"], DALY_VALUE["high"]))
}

point_estimate_ptsd <- function(DW, SEV, PREV, pop_m, pop_f) {
  y <- sum(SEV[c("mild", "moderate", "severe"), "est"] * DW[c("mild", "moderate", "severe"), "est"])
  n <- pop_m * PREV["M", "est"] + pop_f * PREV["F", "est"]
  c(yld_case = y, n_cases = n, yld_total = y * n, rate_100k = y * n / (pop_m + pop_f) * 1e5)
}

# Outputs ----

q3 <- function(x) quantile(x, PROBS, names = FALSE)

summarise_ptsd <- function(sim, hc_2024, cost_2024, axis, condition) {
  mc <- function(label, x) { q <- q3(x)
    data.frame(measure = label, median = q[2], low_95 = q[1], high_95 = q[3]) }
  fixed <- function(label, v) data.frame(measure = label, median = v, low_95 = NA, high_95 = NA)
  out <- rbind(
    mc("DALY per case per year",            sim$yld_case),
    mc("PTSD cases",                        sim$n_cases),
    mc("DALY per year",                     sim$yld_total),
    mc("DALY per 100,000 adults per year",  sim$rate_100k),
    mc("Monetized DALY per case, USD 2024", sim$usd_case),
    fixed("Davis 2022 direct health care cost per case, USD 2024", hc_2024),
    fixed("Davis 2022 total societal cost per case, USD 2024",     cost_2024))
  out$axis <- axis; out$condition <- condition
  out
}

# stability / sensitivity: relative difference (%) of medians and bounds
compare_runs <- function(main, alt, label) {
  m <- c("yld_case", "yld_total", "usd_case")
  do.call(rbind, lapply(m, function(v) {
    a <- q3(main[[v]]); b <- q3(alt[[v]])
    data.frame(check = label, measure = v,
               main_median = a[2], alt_median = b[2],
               diff_median_pct = 100 * (b[2] - a[2]) / a[2],
               diff_low_pct    = 100 * (b[1] - a[1]) / a[1],
               diff_high_pct   = 100 * (b[3] - a[3]) / a[3], row.names = NULL)
  }))
}

plot_ptsd_costs <- function(summ, title, file_stem) {
  lv <- c("Health loss\n(GBD DALY x USD 100-150k)",
          "Direct health care cost\n(Davis et al. 2022)",
          "Total societal cost\n(Davis et al. 2022)")
  get <- function(m) summ[summ$measure == m, ]
  h <- get("Monetized DALY per case, USD 2024")
  pd <- data.frame(estimate = factor(lv, levels = lv),
                   median = c(h$median,
                              get("Davis 2022 direct health care cost per case, USD 2024")$median,
                              get("Davis 2022 total societal cost per case, USD 2024")$median),
                   low  = c(h$low_95, NA, NA), high = c(h$high_95, NA, NA))
  g <- ggplot(pd, aes(estimate, median)) +
    geom_col(width = 0.55, fill = c("grey40", "grey65", "grey80")) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.15, na.rm = TRUE) +
    scale_y_continuous(labels = scales::dollar) +
    labs(x = NULL, y = paste0("USD ", TARGET_YR, " per case per year"), title = title,
         caption = "Error bar: Monte Carlo 95% interval (10,000 draws)") +
    theme_bw()
  ggsave(file.path(FIG_DIR, paste0(file_stem, ".png")), g, width = 7.5, height = 4.5, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0(file_stem, ".pdf")), g, width = 7.5, height = 4.5)
  invisible(g)
}

# Exposures (Axes C and D): attributable fraction and GBD data by age ----

# Relative risk: positive and asymmetric -> LOG-NORMAL.
# A 95% CI of an RR is computed on the log scale, so:
#   meanlog = log(RR), sdlog = (log(high) - log(low)) / (2 * 1.96)
fit_lognormal <- function(est, low, high) {
  est <- unname(est); low <- unname(low); high <- unname(high)   # drop names like "meanlog.est"
  sdlog <- (log(high) - log(low)) / (2 * 1.96)
  q <- qlnorm(PROBS, log(est), sdlog)
  c(meanlog = log(est), sdlog = sdlog, err_pct = rel_err(q, c(low, est, high)))
}
rlnorm_fit <- function(n, fit) rlnorm(n, fit["meanlog"], fit["sdlog"])

# Census margin of error (90% interval) -> 95% interval
moe90_to_ci95 <- function(est, moe90) {
  half <- moe90 / qnorm(.95) * qnorm(.975)
  c(est = est, low = est - half, high = est + half)
}

# Population attributable fraction (Levin): p = prevalence, RR = relative risk
paf <- function(p, rr) p * (rr - 1) / (1 + p * (rr - 1))

# GBD Results Tool (CSV, also in Italian): sums the 5-year age bands into
# age groups. Only the 5-year bands ("25-29", ..., ">=95") are used, not
# the aggregates added by the tool (50-69, 65-74, ...).
# The uncertainty of the sum is approximated by summing the band limits
# (perfect correlation between bands): a conservative interval, still very
# narrow for all-cause mortality estimates.
read_gbd_by_age <- function(file, groups) {
  g <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  g <- g[g$metric_id == 1 & g$measure_id %in% c(1, 2, 4), ]       # 1 = Number; deaths, DALYs, YLL
  lo_age <- suppressWarnings(as.integer(sub("^([0-9]+)-([0-9]+).*$", "\\1", g$age_name)))
  hi_age <- suppressWarnings(as.integer(sub("^([0-9]+)-([0-9]+).*$", "\\2", g$age_name)))
  is95   <- grepl("95", g$age_name) & !grepl("-", g$age_name)
  five   <- (!is.na(lo_age) & !is.na(hi_age) & hi_age - lo_age == 4) | is95
  g$age_start <- ifelse(is95, 95L, lo_age)
  g <- g[five, ]
  measure <- ifelse(g$measure_id == 1, "deaths", ifelse(g$measure_id == 4, "yll",
             ifelse(g$measure_id == 2, "daly", NA)))
  out <- list()
  for (gr in names(groups)) {
    r <- groups[[gr]]
    sel <- g$age_start >= r[1] & g$age_start <= r[2]
    for (m in c("deaths", "yll", "daly")) {
      s <- sel & measure == m
      out[[length(out) + 1]] <- data.frame(group = gr, measure = m,
        est = sum(g$val[s]), low = sum(g$lower[s]), high = sum(g$upper[s]),
        n_age_bands = sum(s))
    }
  }
  do.call(rbind, out)
}

# Population by age group from GBD: number / rate per 100,000
gbd_population_by_age <- function(file, groups) {
  g <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  g <- g[g$measure_id == 1, ]                                  # deaths
  num  <- g[g$metric_id == 1, c("age_name", "val")]
  rate <- g[g$metric_id == 3, c("age_name", "val")]
  m <- merge(num, rate, by = "age_name", suffixes = c("_n", "_r"))
  lo <- suppressWarnings(as.integer(sub("^([0-9]+)-([0-9]+).*$", "\\1", m$age_name)))
  hi <- suppressWarnings(as.integer(sub("^([0-9]+)-([0-9]+).*$", "\\2", m$age_name)))
  is95 <- grepl("95", m$age_name) & !grepl("-", m$age_name)
  keep <- (!is.na(lo) & !is.na(hi) & hi - lo == 4) | is95
  m$age_start <- ifelse(is95, 95L, lo); m <- m[keep, ]
  m$pop <- m$val_n / m$val_r * 1e5
  sapply(groups, function(r) sum(m$pop[m$age_start >= r[1] & m$age_start <= r[2]]))
}
