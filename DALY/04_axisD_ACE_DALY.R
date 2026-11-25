################################################################################
#     AXIS D - BURDEN ATTRIBUTABLE TO ADVERSE CHILDHOOD EXPERIENCES (DALY)
################################################################################
#namefile: DALY/04_axisD_ACE_DALY.R
#creation dt: 13/10/2026
#
# Social cost model, separate from the transcriptomic analysis. The exposure of
# Axis D (GSE70603) is childhood adversity. As for poverty (Axis C), the burden
# ATTRIBUTABLE to it is estimated with the population attributable fraction,
# here for a multi-level exposure (0, 1, 2+ ACEs), following Bellis et al. 2019:
#
#   PAF_k = [p1 (RR1_k - 1) + p2 (RR2_k - 1)] / [1 + p1 (RR1_k - 1) + p2 (RR2_k - 1)]
#   (same as the Bellis formula: p0 + p1 RR1 + p2 RR2 = 1 + p1(RR1-1) + p2(RR2-1))
#   attributable DALYs_k = PAF_k * DALYs_k (GBD 2023, US, adults aged 25+)
#   total = sum over the 6 causes (GBD causes are mutually exclusive)
#
# The 6 causes of Bellis: anxiety, depression, neoplasms, cardiovascular
# diseases, diabetes, chronic respiratory diseases. The 4 risk factors of Bellis
# (alcohol, drugs, smoking, obesity) are NOT included: they overlap with the
# causes (e.g. smoking causes cancer) and would need double counting to be
# removed by hand, so the estimate is conservative.
# Unlike Axis C, full DALYs (YLL + YLD) are obtained, because GBD causes have
# both components; the RRs of Bellis do not depend on age, so the total for
# ages 25+ is enough.
#
# Inputs:
#   - pooled RRs, North America, 1 ACE and 2+ ACEs vs 0 ACEs, by cause
#     [Bellis et al. 2019, Lancet Public Health 4(10):e517-e528, Table 2]
#   - pooled ACE prevalence, North America (9 studies): 1 ACE 23.4% (22.0-24.8),
#     2+ ACEs 35.0% (31.6-38.4) [Bellis 2019, Results]
#   - DALYs, YLL and YLD by cause, US, 2023, aged 25+
#     [GBD 2023 Results; data/source/GBD2023_CAUSE_6causes_AXIS_D/]
#
# External checks:
#   (a) DALYs attributable to "child maltreatment" in GBD 2023 (childhood
#       sexual abuse and bullying), US, aged 25+
#       [data/source/GBD2023_RISK_child_maltreatment_AXIS_D/]: a narrower
#       measure than ACEs, so a lower value is expected
#   (b) Bellis 2019: 8.40 million DALYs attributable to ACEs for the 6 causes
#       in North America (US and Canada), 2017, aged 15+ (sum of Table 3)
#
# Uncertainty: Monte Carlo, 10,000 draws, seed 18: log-normal RRs, beta
# prevalences (quantile matching), log-normal GBD DALYs, uniform value of a
# DALY. Stability with seed 2026.
#
# Outputs:
#   results/D_ACE_DALY.csv, results/D_ACE_DALY_by_cause.csv,
#   results/D_distribution_fit.csv, results/D_sensitivity.csv,
#   figures/D_ACE_DALY_by_cause.png/.pdf
#
# Run from the project root:
#   source("DALY/04_axisD_ACE_DALY.R", echo = TRUE)
################################################################################

source("DALY/00_daly_functions.R")

AXIS      <- "D"
CONDITION <- "Adverse childhood experiences (ACEs), US adults 25+"
N_SIM     <- 10000

# Bellis 2019, Table 2, North America: RR (95% CI) for 1 ACE and 2+ ACEs
RR1 <- rbind(anxiety     = c(est = 1.08, low = 0.66, high = 1.77),
             depression  = c(est = 1.34, low = 1.15, high = 1.56),
             cancer      = c(est = 1.10, low = 0.95, high = 1.28),
             cvd         = c(est = 1.16, low = 1.04, high = 1.29),
             diabetes    = c(est = 1.11, low = 1.02, high = 1.22),
             respiratory = c(est = 1.29, low = 1.17, high = 1.41))
RR2 <- rbind(anxiety     = c(est = 2.25, low = 1.43, high = 3.56),
             depression  = c(est = 2.69, low = 2.17, high = 3.33),
             cancer      = c(est = 1.25, low = 1.10, high = 1.43),
             cvd         = c(est = 1.60, low = 1.41, high = 1.81),
             diabetes    = c(est = 1.15, low = 0.97, high = 1.37),
             respiratory = c(est = 1.90, low = 1.58, high = 2.29))
# Bellis 2019, pooled prevalence, North America
PREV <- rbind(ace1  = c(est = .234, low = .220, high = .248),
              ace2p = c(est = .350, low = .316, high = .384))

# GBD cause names (download in Italian) -> short names
CAUSE_MAP <- c("Disturbi d'ansia" = "anxiety", "Disturbi depressivi" = "depression",
               "Neoplasie" = "cancer", "Malattie cardiovascolari" = "cvd",
               "Diabete mellito" = "diabetes", "Malattie respiratorie croniche" = "respiratory")
CAUSE_LABEL <- c(anxiety = "Anxiety disorders", depression = "Depressive disorders",
                 cancer = "Neoplasms", cvd = "Cardiovascular diseases",
                 diabetes = "Diabetes mellitus", respiratory = "Chronic respiratory diseases")
BELLIS_NAM_DALY <- (622 + 798 + 1641 + 3321 + 313 + 1707) * 1e3   # Table 3, thousands

one_csv <- function(dir) {
  f <- list.files(file.path(SRC_DIR, dir), "\\.csv$", full.names = TRUE)
  stopifnot(length(f) == 1); f
}

# GBD inputs ----
g <- read.csv(one_csv("GBD2023_CAUSE_6causes_AXIS_D"), stringsAsFactors = FALSE)
g$cause <- CAUSE_MAP[g$cause_name]
stopifnot(!any(is.na(g$cause)), all(g$age_name == unique(g$age_name)[1]))
cat("Ages:", unique(g$age_name), "| year:", unique(g$year), "\n")
get_m <- function(m_id) {                           # 2 = DALY, 3 = YLD, 4 = YLL
  x <- g[g$measure_id == m_id & g$metric_id == 1, ]
  out <- t(sapply(names(CAUSE_LABEL), function(k) {
    r <- x[x$cause == k, ]
    if (nrow(r) == 0) return(c(est = 0, low = 0, high = 0))    # e.g. YLL of anxiety / depression
    c(est = r$val, low = r$lower, high = r$upper)
  }))
  out
}
DALY <- get_m(2); YLD <- get_m(3); YLL <- get_m(4)
print(round(DALY))
d_rate <- g[g$measure_id == 2 & g$metric_id == 3, ]
POP <- median(d_rate$val[match(names(CAUSE_LABEL), d_rate$cause)] ^ -1 *
              DALY[, "est"] * 1e5)                              # number / rate
cat("US population aged 25+ (GBD 2023):", format(round(POP), big.mark = ","), "\n")

rk <- read.csv(one_csv("GBD2023_RISK_child_maltreatment_AXIS_D"), stringsAsFactors = FALSE)
rk <- rk[rk$metric_id == 1 & rk$measure_id == 2, c("rei_name", "val", "lower", "upper")]
print(rk)
GBD_MALTREAT <- rk[rk$rei_name == "Maltrattamento di minori", ]

# Distributions: published vs reproduced ----
fit_rr1 <- lapply(rownames(RR1), function(k) fit_lognormal(RR1[k, 1], RR1[k, 2], RR1[k, 3])); names(fit_rr1) <- rownames(RR1)
fit_rr2 <- lapply(rownames(RR2), function(k) fit_lognormal(RR2[k, 1], RR2[k, 2], RR2[k, 3])); names(fit_rr2) <- rownames(RR2)
fit_p   <- lapply(rownames(PREV), function(k) fit_beta(PREV[k, 1], PREV[k, 2], PREV[k, 3]));  names(fit_p)   <- rownames(PREV)
fit_dl  <- lapply(rownames(DALY), function(k) fit_lognormal(DALY[k, 1], DALY[k, 2], DALY[k, 3])); names(fit_dl) <- rownames(DALY)
fits <- rbind(
  data.frame(group = "RR 1 ACE (lognormal)",  parameter = rownames(RR1), RR1, err_pct = sapply(fit_rr1, `[`, "err_pct"), row.names = NULL),
  data.frame(group = "RR >=2 ACE (lognormal)", parameter = rownames(RR2), RR2, err_pct = sapply(fit_rr2, `[`, "err_pct"), row.names = NULL),
  data.frame(group = "ACE prevalence (beta)", parameter = rownames(PREV), PREV, err_pct = sapply(fit_p, `[`, "err_pct"), row.names = NULL),
  data.frame(group = "GBD DALY (lognormal)", parameter = rownames(DALY), DALY, err_pct = sapply(fit_dl, `[`, "err_pct"), row.names = NULL))
print(fits, digits = 3)
write.csv(fits, file.path(RES_DIR, paste0(AXIS, "_distribution_fit.csv")), row.names = FALSE)

# Point estimate ----
paf_multi <- function(p1, p2, rr1, rr2) {
  x <- p1 * (rr1 - 1) + p2 * (rr2 - 1); x / (1 + x)
}
pt <- data.frame(cause = rownames(DALY),
                 paf = paf_multi(PREV["ace1", "est"], PREV["ace2p", "est"], RR1[, "est"], RR2[, "est"]),
                 daly = DALY[, "est"])
pt$daly_attr <- pt$paf * pt$daly
print(pt, digits = 3)
cat("Total attributable DALYs (point estimate):", format(round(sum(pt$daly_attr)), big.mark = ","), "\n")

# Monte Carlo ----
run_ace <- function(seed) {
  set.seed(seed)
  p1 <- rbeta_fit(N_SIM, fit_p[["ace1"]]); p2 <- rbeta_fit(N_SIM, fit_p[["ace2p"]])
  by_cause <- lapply(rownames(DALY), function(k) {
    f <- paf_multi(p1, p2, rlnorm_fit(N_SIM, fit_rr1[[k]]), rlnorm_fit(N_SIM, fit_rr2[[k]]))
    d <- rlnorm_fit(N_SIM, fit_dl[[k]])
    share_yll <- YLL[k, "est"] / DALY[k, "est"]          # YLL share of the cause (GBD estimate)
    data.frame(paf = f, daly_attr = f * d, yll_attr = f * d * share_yll,
               yld_attr = f * d * (1 - share_yll))
  })
  names(by_cause) <- rownames(DALY)
  tot <- data.frame(daly_attr = Reduce(`+`, lapply(by_cause, `[[`, "daly_attr")),
                    yll_attr  = Reduce(`+`, lapply(by_cause, `[[`, "yll_attr")),
                    yld_attr  = Reduce(`+`, lapply(by_cause, `[[`, "yld_attr")),
                    exposed   = (p1 + p2) * POP)
  tot$daly_per_exposed <- tot$daly_attr / tot$exposed
  tot$rate_100k        <- tot$daly_attr / POP * 1e5
  tot$usd_per_exposed  <- tot$daly_per_exposed * runif(N_SIM, DALY_VALUE["low"], DALY_VALUE["high"])
  tot$usd_total        <- tot$daly_attr * runif(N_SIM, DALY_VALUE["low"], DALY_VALUE["high"])
  list(by_cause = by_cause, total = tot)
}
sim    <- run_ace(18)
sim_s2 <- run_ace(2026)

mc <- function(label, x) { q <- q3(x)
  data.frame(measure = label, median = q[2], low_95 = q[1], high_95 = q[3]) }
summ <- rbind(
  mc("Adults with >=1 ACE, 25+",                      sim$total$exposed),
  mc("Attributable DALY per year, 6 causes",          sim$total$daly_attr),
  mc("  of which YLL",                                sim$total$yll_attr),
  mc("  of which YLD",                                sim$total$yld_attr),
  mc("Attributable DALY per 100,000 adults 25+ per year", sim$total$rate_100k),
  mc("Attributable DALY per person with >=1 ACE per year", sim$total$daly_per_exposed),
  mc("Monetized DALY per person with >=1 ACE, USD 2024",   sim$total$usd_per_exposed),
  mc("Monetized attributable DALY per year, USD 2024",     sim$total$usd_total),
  data.frame(measure = "Check (a): GBD 2023 DALY attributable to child maltreatment",
             median = GBD_MALTREAT$val, low_95 = GBD_MALTREAT$lower, high_95 = GBD_MALTREAT$upper),
  data.frame(measure = "Check (b): Bellis 2019 DALY attributable to ACEs, 6 causes, North America 2017",
             median = BELLIS_NAM_DALY, low_95 = NA, high_95 = NA))
summ$axis <- AXIS; summ$condition <- CONDITION
print(summ, digits = 4)
write.csv(summ, file.path(RES_DIR, paste0(AXIS, "_ACE_DALY.csv")), row.names = FALSE)

bc <- do.call(rbind, lapply(names(sim$by_cause), function(k) {
  a <- q3(sim$by_cause[[k]]$paf); b <- q3(sim$by_cause[[k]]$daly_attr)
  data.frame(cause = CAUSE_LABEL[k], paf_median = a[2], paf_low = a[1], paf_high = a[3],
             daly_total = DALY[k, "est"], daly_attr_median = b[2], daly_attr_low = b[1], daly_attr_high = b[3],
             row.names = NULL)
}))
bc <- bc[order(-bc$daly_attr_median), ]
print(bc, digits = 3)
write.csv(bc, file.path(RES_DIR, paste0(AXIS, "_ACE_DALY_by_cause.csv")), row.names = FALSE)

sens <- do.call(rbind, lapply(c("daly_attr", "daly_per_exposed", "usd_per_exposed"), function(v) {
  a <- q3(sim$total[[v]]); b <- q3(sim_s2$total[[v]])
  data.frame(check = "stability: seed 2026 vs 18", measure = v,
             main_median = a[2], alt_median = b[2],
             diff_median_pct = 100 * (b[2] - a[2]) / a[2],
             diff_low_pct = 100 * (b[1] - a[1]) / a[1],
             diff_high_pct = 100 * (b[3] - a[3]) / a[3])
}))
print(sens, digits = 3)
write.csv(sens, file.path(RES_DIR, paste0(AXIS, "_sensitivity.csv")), row.names = FALSE)

# Plot: attributable DALYs by cause ----
bc$cause <- factor(bc$cause, levels = rev(bc$cause))
gp <- ggplot(bc, aes(cause, daly_attr_median / 1e6)) +
  geom_col(width = 0.6, fill = "grey55") +
  geom_errorbar(aes(ymin = daly_attr_low / 1e6, ymax = daly_attr_high / 1e6), width = 0.2) +
  coord_flip() +
  labs(x = NULL, y = "Attributable DALYs per year (millions)",
       title = "Axis D - DALYs attributable to ACEs, US adults 25+, 2023",
       caption = "Error bar: Monte Carlo 95% interval (10,000 draws)") +
  theme_bw()
ggsave(file.path(FIG_DIR, "D_ACE_DALY_by_cause.png"), gp, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "D_ACE_DALY_by_cause.pdf"), gp, width = 7, height = 4.5)

# Check with Peterson et al. 2018 (child maltreatment) ----
# Peterson, Florence & Klevens 2018, Child Abuse Negl 86:178-183. Two
# comparisons, both only indicative:
#   (i) health units: loss of quality of life in adults reporting childhood
#       maltreatment = 0.028 (0.022-0.034) SF-6D utility per year of adult
#       life, i.e. QALYs per year [Corso et al. 2008, Am J Public Health
#       98(6):1094-1100; used by Peterson 2018]. It is compared with our
#       attributable YLD per person with 1+ ACE per year (morbidity only, like
#       the QALY of Corso). The exposures differ: maltreatment (abuse and
#       neglect) vs 1+ ACE (which includes household dysfunction).
#   (ii) monetary units: lifetime cost per (non-fatal) victim $830,928 and
#       annual US burden $428 billion (substantiated cases, 2015), in USD 2015
#       -> USD 2024 with CPI-U. Peterson values health with the VSL ($16.6
#       million) and discounts at 3%: not comparable with our value per DALY
#       ($100-150k, annual), so it is reported as context.
PET_QALY_YR   <- 0.028
PET_VICTIM    <- 830928
PET_BURDEN    <- 428e9
PET_YR        <- 2015
f15 <- inflation_factor(read_cpi(), PET_YR)
yld_pe <- q3(sim$total$yld_attr / sim$total$exposed)
pet <- data.frame(
  measure = c("Attributable YLD per person with >=1 ACE per year (this study)",
              "Corso 2008 / Peterson 2018: QALY loss per year of adulthood, maltreatment",
              "Peterson 2018: lifetime cost per nonfatal victim, USD 2024",
              "Peterson 2018: annual US economic burden (substantiated), USD 2024"),
  median  = c(yld_pe[2], PET_QALY_YR, PET_VICTIM * f15, PET_BURDEN * f15),
  low_95  = c(yld_pe[1], 0.022, NA, NA),
  high_95 = c(yld_pe[3], 0.034, NA, NA))
cat("CPI-U factor 2015 -> 2024:", round(f15, 4), "\n")
print(pet, digits = 4)
write.csv(pet, file.path(RES_DIR, paste0(AXIS, "_check_Peterson.csv")), row.names = FALSE)
