################################################################################
#                AXIS C - BURDEN ATTRIBUTABLE TO POVERTY (DALY)
################################################################################
#namefile: DALY/03_axisC_poverty_DALY.R
#creation dt: 13/10/2026
#
# Social cost model, separate from the transcriptomic analysis. The exposure of
# Axis C (GSE221766, HANDLS) is poverty, which is not a disease and has no
# disability weights. The burden ATTRIBUTABLE to poverty is therefore estimated
# with the population attributable fraction (PAF), as in Galea et al. 2011, but
# on years of life lost (YLL) instead of deaths only:
#
#   PAF_g = p_g (RR_g - 1) / (1 + p_g (RR_g - 1))          for age group g
#   attributable YLL_g = PAF_g * all-cause YLL_g (GBD 2023, US)
#   ages 25-64 and 65+ (different RRs by age in Galea 2011)
#
# Only mortality (YLL) is counted: the contribution of poverty to morbidity
# (YLD) is not included, so the estimate is conservative. For this reason Axis C
# is expressed in attributable YLL, which are a part of the DALYs.
#
# Inputs:
#   - RR of all-cause mortality, poor vs non-poor (meta-analysis, individual
#     poverty around the federal threshold): ages 25-64 1.75 (1.51-2.04); 65+
#     1.40 (1.37-1.43)  [Galea et al. 2011, AJPH 101(8):1456-1465, Table 3]
#   - official poverty rate (100% threshold), 2023: ages 18-64 10.0% (90% MOE
#     0.3), 65+ 9.7% (0.4)  [Shrider & Bijou 2025, Census P60-287, Table A-1].
#     The 100% threshold (not the 125% of HANDLS) matches the definition behind
#     the RRs of Galea.
#   - all-cause deaths, YLL and DALYs, US, 2023, 5-year age bands
#     [GBD 2023 Results, IHME 2024; data/source/GBD2023_CAUSE_allcauses_AXIS_C/]
#
# External check: attributable deaths vs Galea 2011 (133,250 in 2000).
#
# Uncertainty: Monte Carlo, 10,000 draws, seed 18: log-normal RR (CI symmetric
# on the log scale), beta prevalences (quantile matching), log-normal GBD YLL
# and deaths (GBD interval), uniform value of a DALY. Stability with seed 2026.
#
# Outputs:
#   results/C_poverty_DALY.csv, results/C_distribution_fit.csv,
#   results/C_sensitivity.csv, figures/C_poverty_YLL_by_age.png/.pdf
#
# Run from the project root:
#   source("DALY/03_axisC_poverty_DALY.R", echo = TRUE)
################################################################################

source("DALY/00_daly_functions.R")

AXIS      <- "C"
CONDITION <- "Poverty (official threshold), US adults 25+"
GROUPS    <- list(`25-64` = c(25, 64), `65+` = c(65, 95))
N_SIM     <- 10000

RR <- rbind(`25-64` = c(est = 1.75, low = 1.51, high = 2.04),   # Galea 2011
            `65+`   = c(est = 1.40, low = 1.37, high = 1.43))
PREV <- rbind(`25-64` = moe90_to_ci95(.100, .003),              # Census 2023, 18-64
              `65+`   = moe90_to_ci95(.097, .004))              # Census 2023, 65+
GALEA_DEATHS <- 133250                                           # attributable deaths, 2000

# folder of the GBD 2023 all-cause download for Axis C
gbd_file <- list.files(file.path(SRC_DIR, "GBD2023_CAUSE_allcauses_AXIS_C"), "\\.csv$",
                       full.names = TRUE)
stopifnot(length(gbd_file) == 1)
cat("File GBD:", gbd_file, "\n")

# Inputs ----
gbd <- read_gbd_by_age(gbd_file, GROUPS)
print(gbd)
stopifnot(all(gbd$n_age_bands[gbd$group == "25-64"] == 8),
          all(gbd$n_age_bands[gbd$group == "65+"]   == 7))
pop <- gbd_population_by_age(gbd_file, GROUPS)
cat("GBD 2023 population:", paste(names(pop), format(round(pop), big.mark = ","), collapse = " | "), "\n")
print(RR); print(PREV)

# Distributions: published vs reproduced ----
fit_rr <- lapply(rownames(RR), function(g) fit_lognormal(RR[g, "est"], RR[g, "low"], RR[g, "high"]))
names(fit_rr) <- rownames(RR)
fit_p  <- lapply(rownames(PREV), function(g) fit_beta(PREV[g, "est"], PREV[g, "low"], PREV[g, "high"]))
names(fit_p) <- rownames(PREV)
get_g <- function(g, m) gbd[gbd$group == g & gbd$measure == m, ]
fit_gbd <- list()
for (g in names(GROUPS)) for (m in c("deaths", "yll")) {
  r <- get_g(g, m); fit_gbd[[paste(g, m)]] <- fit_lognormal(r$est, r$low, r$high)
}
fits <- rbind(
  data.frame(group = "relative risk (lognormal)", parameter = rownames(RR), RR,
             err_pct = sapply(fit_rr, `[`, "err_pct"), row.names = NULL),
  data.frame(group = "poverty prevalence (beta)", parameter = rownames(PREV), PREV,
             err_pct = sapply(fit_p, `[`, "err_pct"), row.names = NULL))
print(fits, digits = 3)
write.csv(fits, file.path(RES_DIR, paste0(AXIS, "_distribution_fit.csv")), row.names = FALSE)

# Point estimate ----
pt <- do.call(rbind, lapply(names(GROUPS), function(g) {
  f <- paf(PREV[g, "est"], RR[g, "est"])
  data.frame(group = g, paf = f,
             deaths_attr = f * get_g(g, "deaths")$est,
             yll_attr    = f * get_g(g, "yll")$est,
             exposed     = PREV[g, "est"] * pop[g])
}))
pt$yll_per_exposed <- pt$yll_attr / pt$exposed
print(pt, digits = 4)
cat("Total 25+: attributable deaths", round(sum(pt$deaths_attr)),
    "| attributable YLL", round(sum(pt$yll_attr)),
    "| YLL per person in poverty", round(sum(pt$yll_attr) / sum(pt$exposed), 4), "\n")

# Monte Carlo ----
run_poverty <- function(seed) {
  set.seed(seed)
  res <- lapply(names(GROUPS), function(g) {
    p   <- rbeta_fit(N_SIM, fit_p[[g]])
    rr  <- rlnorm_fit(N_SIM, fit_rr[[g]])
    f   <- paf(p, rr)
    dth <- rlnorm_fit(N_SIM, fit_gbd[[paste(g, "deaths")]])
    yll <- rlnorm_fit(N_SIM, fit_gbd[[paste(g, "yll")]])
    data.frame(paf = f, deaths_attr = f * dth, yll_attr = f * yll, exposed = p * pop[g])
  })
  names(res) <- names(GROUPS)
  tot <- data.frame(
    deaths_attr = res[[1]]$deaths_attr + res[[2]]$deaths_attr,
    yll_attr    = res[[1]]$yll_attr + res[[2]]$yll_attr,
    exposed     = res[[1]]$exposed + res[[2]]$exposed)
  tot$yll_per_exposed <- tot$yll_attr / tot$exposed
  tot$rate_100k       <- tot$yll_attr / sum(pop) * 1e5
  tot$usd_per_exposed <- tot$yll_per_exposed * runif(N_SIM, DALY_VALUE["low"], DALY_VALUE["high"])
  list(by_age = res, total = tot)
}
sim    <- run_poverty(18)
sim_s2 <- run_poverty(2026)

mc <- function(label, x) { q <- q3(x)
  data.frame(measure = label, median = q[2], low_95 = q[1], high_95 = q[3]) }
summ <- rbind(
  mc("PAF of all-cause mortality, 25-64", sim$by_age[["25-64"]]$paf),
  mc("PAF of all-cause mortality, 65+",   sim$by_age[["65+"]]$paf),
  mc("Adults in poverty, 25+",            sim$total$exposed),
  mc("Attributable deaths per year, 25+", sim$total$deaths_attr),
  mc("Attributable YLL per year, 25-64",  sim$by_age[["25-64"]]$yll_attr),
  mc("Attributable YLL per year, 65+",    sim$by_age[["65+"]]$yll_attr),
  mc("Attributable YLL per year, 25+",    sim$total$yll_attr),
  mc("Attributable YLL per 100,000 adults 25+ per year", sim$total$rate_100k),
  mc("Attributable YLL per person in poverty per year",  sim$total$yll_per_exposed),
  mc("Monetized YLL per person in poverty, USD 2024",    sim$total$usd_per_exposed))
summ$axis <- AXIS; summ$condition <- CONDITION
print(summ, digits = 4)
write.csv(summ, file.path(RES_DIR, paste0(AXIS, "_poverty_DALY.csv")), row.names = FALSE)

cat("Check: attributable deaths", round(q3(sim$total$deaths_attr)[2]),
    "vs Galea 2011 (year 2000):", GALEA_DEATHS, "\n")

sens <- do.call(rbind, lapply(c("deaths_attr", "yll_attr", "yll_per_exposed", "usd_per_exposed"), function(v) {
  a <- q3(sim$total[[v]]); b <- q3(sim_s2$total[[v]])
  data.frame(check = "stability: seed 2026 vs 18", measure = v,
             main_median = a[2], alt_median = b[2],
             diff_median_pct = 100 * (b[2] - a[2]) / a[2],
             diff_low_pct = 100 * (b[1] - a[1]) / a[1],
             diff_high_pct = 100 * (b[3] - a[3]) / a[3])
}))
print(sens, digits = 3)
write.csv(sens, file.path(RES_DIR, paste0(AXIS, "_sensitivity.csv")), row.names = FALSE)

# Plot: attributable YLL by age ----
pd <- data.frame(age = factor(names(GROUPS), levels = names(GROUPS)),
                 median = sapply(names(GROUPS), function(g) q3(sim$by_age[[g]]$yll_attr)[2]),
                 low    = sapply(names(GROUPS), function(g) q3(sim$by_age[[g]]$yll_attr)[1]),
                 high   = sapply(names(GROUPS), function(g) q3(sim$by_age[[g]]$yll_attr)[3]))
g <- ggplot(pd, aes(age, median / 1e6)) +
  geom_col(width = 0.5, fill = "grey55") +
  geom_errorbar(aes(ymin = low / 1e6, ymax = high / 1e6), width = 0.15) +
  labs(x = "Age group (years)", y = "Attributable YLL per year (millions)",
       title = "Axis C - Years of life lost attributable to poverty, US 2023",
       caption = "Error bar: Monte Carlo 95% interval (10,000 draws)") +
  theme_bw()
ggsave(file.path(FIG_DIR, "C_poverty_YLL_by_age.png"), g, width = 6.5, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "C_poverty_YLL_by_age.pdf"), g, width = 6.5, height = 4.5)

# Sensitivity: current vs cumulative poverty (Brady et al. 2023) ----
# The RRs of Galea 2011 measure poverty at a single point in time. Brady, Kohler
# & Zheng (preprint, 22/01/2023; PSID 1997-2019, N = 18,995, aged 15+, Cox
# models with pre-exposure controls, Model 3):
#   - relative poverty: equivalized disposable income < 50% of the median
#   - current poverty:    HR 1.42 (1.26-1.60)
#   - cumulative poverty (always vs never poor in the last 10 years):
#                         HR 1.71 (1.45-2.02)
#   - relative poverty rate, aged 15+, 2019: 16.4% (LIS), used by Brady for
#     BOTH PAFs (6.5% and 10.5%); no interval, so it is fixed
#   - check: 183,003 (current) and 295,431 (cumulative) deaths, US 2019
# A single HR is applied to all adults aged 25+ (Brady gives no HR by age),
# with the same GBD 2023 YLL and deaths. Differences from the main estimate:
# definition (relative, 50% of the median, vs official threshold), prevalence
# (16.4% vs about 10%), HR not by age. The cumulative estimate uses the
# prevalence of current poverty (as Brady does): it is an upper bound, since
# the always poor are fewer than the poor in a given year.
BRADY_HR <- rbind(current    = c(est = 1.42, low = 1.26, high = 1.60),
                  cumulative = c(est = 1.71, low = 1.45, high = 2.02))
BRADY_P      <- 0.164
BRADY_DEATHS <- c(current = 183003, cumulative = 295431)
pop25 <- sum(pop)

run_brady <- function(seed, hr) {
  set.seed(seed)
  fit <- fit_lognormal(unname(hr["est"]), unname(hr["low"]), unname(hr["high"]))
  rr  <- rlnorm_fit(N_SIM, fit)
  f   <- paf(BRADY_P, rr)
  dth <- Reduce(`+`, lapply(names(GROUPS), function(g) rlnorm_fit(N_SIM, fit_gbd[[paste(g, "deaths")]])))
  yll <- Reduce(`+`, lapply(names(GROUPS), function(g) rlnorm_fit(N_SIM, fit_gbd[[paste(g, "yll")]])))
  data.frame(paf = f, deaths_attr = f * dth, yll_attr = f * yll,
             yll_per_exposed = f * yll / (BRADY_P * pop25),
             rate_100k = f * yll / pop25 * 1e5)
}
brady <- do.call(rbind, lapply(rownames(BRADY_HR), function(s) {
  b <- run_brady(18, BRADY_HR[s, ])
  do.call(rbind, lapply(c("paf", "deaths_attr", "yll_attr", "rate_100k", "yll_per_exposed"), function(v) {
    q <- q3(b[[v]]); data.frame(scenario = paste("Brady 2023,", s, "poverty"), measure = v,
                                median = q[2], low_95 = q[1], high_95 = q[3])
  }))
}))
main <- do.call(rbind, lapply(c("deaths_attr", "yll_attr", "rate_100k", "yll_per_exposed"), function(v) {
  q <- q3(sim$total[[v]]); data.frame(scenario = "Main: Galea 2011, official poverty", measure = v,
                                      median = q[2], low_95 = q[1], high_95 = q[3])
}))
brady <- rbind(main, brady)
rownames(brady) <- NULL
print(brady, digits = 4)
cat("Brady check: attributable deaths 2019 (aged 15+), current", BRADY_DEATHS["current"],
    "| cumulative", BRADY_DEATHS["cumulative"], "\n")
cat("YLL ratio, cumulative / current:",
    round(brady$median[brady$scenario == "Brady 2023, cumulative poverty" & brady$measure == "yll_attr"] /
          brady$median[brady$scenario == "Brady 2023, current poverty" & brady$measure == "yll_attr"], 2), "\n")
write.csv(brady, file.path(RES_DIR, paste0(AXIS, "_sensitivity_Brady.csv")), row.names = FALSE)
