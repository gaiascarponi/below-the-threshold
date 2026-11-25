################################################################################
#                AXIS A - BURDEN OF PTSD (DALY) - US CIVILIANS
################################################################################
#namefile: DALY/01_axisA_PTSD_DALY.R
#creation dt: 13/10/2026
#
# Social cost model, separate from the transcriptomic analysis: it does not use
# the expression data. It estimates the burden of PTSD (the clinical outcome of
# Axis A, GSE97356) in the US adult civilian population, in DALYs, with the
# Global Burden of Disease method.
#
#   DALY = YLL + YLD;  YLL = 0 (GBD assigns no direct mortality to anxiety
#   disorders);  YLD per case = sum_s p_s * DW_s
#   PTSD is not a separate GBD cause: it is part of "anxiety disorders"
#   (ICD-10 F43.1, DSM-IV-TR 309.81; Yang et al. 2021), so the disability
#   weights and severity distribution of that cause are used.
#
# Inputs:
#   - GBD 2019 disability weights (IHME 2020), file in data/source
#   - GBD 2019 severity distribution of anxiety disorders (Yang et al. 2021,
#     Epidemiol Psychiatr Sci 30:e36, supplementary) -> 00_daly_functions.R
#   - 12-month PTSD prevalence, civilians (DSM-5, NESARC-III 2012-13): M 2.6%
#     (2.2-2.9), F 6.0% (5.5-6.6) [Lehavot et al. 2018, Am J Prev Med 54:e1-e9]
#   - US adult civilian populations 2018 and costs per case (USD 2018) from
#     Davis et al. 2022, J Clin Psychiatry 83(3):21m14116: 2,607,131 M and
#     7,558,833 F with PTSD; total societal cost $18,640; direct health care
#     costs $6,495 (Table 1)
#   - CPI-U CUUR0000SA0 (BLS), file in data/source
#
# Uncertainty: Monte Carlo, 10,000 draws, seed 18 (00_daly_functions.R): beta
#   (quantile matching) for weights, severity and prevalences; uniform for the
#   value of a DALY. Checks: (i) stability with seed 2026; (ii) sensitivity
#   with a Dirichlet for severity.
#
# Outputs:
#   results/A_PTSD_DALY.csv               main estimates
#   results/A_distribution_fit.csv        published vs reproduced (beta / normal)
#   results/A_sensitivity.csv             stability and sensitivity
#   figures/A_PTSD_DALY_vs_cost.png/.pdf
#
# Run from the project root:
#   source("DALY/01_axisA_PTSD_DALY.R", echo = TRUE)
################################################################################

source("DALY/00_daly_functions.R")

AXIS      <- "A"
CONDITION <- "PTSD (US civilian adults)"
POP_M     <- 2607131 / 0.026          # US adult civilians 2018, men (Davis 2022)
POP_F     <- 7558833 / 0.060          # US adult civilians 2018, women (Davis 2022)
PREV      <- rbind(M = c(est = .026, low = .022, high = .029),     # Lehavot 2018
                   F = c(est = .060, low = .055, high = .066))
COST_TOT  <- 18640                    # USD 2018 per case, total societal cost
COST_HC   <- 6495                     # USD 2018 per case, direct health care costs
COST_YR   <- 2018

# Inputs ----
DW <- read_disability_weights(c(mild     = "Mild anxiety disorders",
                                moderate = "Moderate anxiety disorders",
                                severe   = "Severe anxiety disorders"))
print(DW)                                     # expected 0.030 / 0.133 / 0.523
SEV <- SEVERITY_ANXIETY
print(SEV)

cpi <- read_cpi()
f   <- inflation_factor(cpi, COST_YR)
cat("CPI-U annual mean", COST_YR, "=", round(cpi$annual[cpi$year == COST_YR], 3),
    "|", TARGET_YR, "=", round(cpi$annual[cpi$year == TARGET_YR], 3),
    "| factor =", round(f, 4), "\n")

# Distributions: published vs reproduced ----
fits <- rbind(fit_table(DW, "disability weight"),
              fit_table(SEV, "severity proportion"),
              fit_table(PREV, "prevalence"))
print(fits[, c("group", "parameter", "published_low", "published_est", "published_high",
               "beta_q025", "beta_q50", "beta_q975", "beta_err_pct", "truncnorm_err_pct")],
      digits = 3)
dir_fit <- fit_dirichlet(SEV)
cat("Dirichlet (sensitivity): K =", round(dir_fit$K, 1),
    "| error % by category:", paste(names(dir_fit$err_pct), round(dir_fit$err_pct, 1), collapse = ", "), "\n")
write.csv(fits, file.path(RES_DIR, paste0(AXIS, "_distribution_fit.csv")), row.names = FALSE)

# Point estimate ----
pt <- point_estimate_ptsd(DW, SEV, PREV, POP_M, POP_F)
print(round(pt, 4))

# Monte Carlo (main analysis) and checks ----
sim     <- run_ptsd_daly(DW, SEV, PREV, POP_M, POP_F, seed = 18)
sim_s2  <- run_ptsd_daly(DW, SEV, PREV, POP_M, POP_F, seed = 2026)
sim_dir <- run_ptsd_daly(DW, SEV, PREV, POP_M, POP_F, seed = 18, severity_method = "dirichlet")

summ <- summarise_ptsd(sim, COST_HC * f, COST_TOT * f, AXIS, CONDITION)
print(summ, digits = 4)
write.csv(summ, file.path(RES_DIR, paste0(AXIS, "_PTSD_DALY.csv")), row.names = FALSE)

sens <- rbind(compare_runs(sim, sim_s2,  "stability: seed 2026 vs 18"),
              compare_runs(sim, sim_dir, "sensitivity: Dirichlet severity"))
print(sens, digits = 3)
write.csv(sens, file.path(RES_DIR, paste0(AXIS, "_sensitivity.csv")), row.names = FALSE)

# Plot ----
plot_ptsd_costs(summ, "Axis A - PTSD (US civilians): annual burden per case",
                paste0(AXIS, "_PTSD_DALY_vs_cost"))

# Check with GBD 2023 and update to 2023 (adults aged 25+) ----
# (i) Independent check of the DALY per case: for the GBD cause "anxiety
#     disorders" (which includes PTSD), YLD / prevalent cases, US 2023, aged
#     25+, by sex [GBD 2023; data/source/GBD2023_CAUSE_anxiety_AXIS_A/]. This is
#     the DALY per case computed by GBD itself, with GBD 2023 weights and severity.
# (ii) Update to 2023 on the same basis as Axes C and D (US, aged 25+, GBD 2023
#     population): cases = pop_M * prev_M + pop_F * prev_F, with the Lehavot 2018
#     prevalences. The GBD population is the total one (it includes the
#     military, about 1% of adults): it is used for the comparison across axes,
#     while the main estimate above stays on the civilians of Davis 2022,
#     consistent with the costs.
anx <- read.csv(list.files(file.path(SRC_DIR, "GBD2023_CAUSE_anxiety_AXIS_A"), "\\.csv$",
                           full.names = TRUE), stringsAsFactors = FALSE)
stopifnot(length(unique(anx$age_name)) == 1, unique(anx$year) == 2023)
val <- function(measure_id, metric_id, sex_id, col = "val")
  anx[anx$measure_id == measure_id & anx$metric_id == metric_id & anx$sex_id == sex_id, col]
# 3 = YLD, 5 = prevalence; metric 1 = number, 3 = rate per 100,000; sex 1 = M, 2 = F
gbd_check <- data.frame(
  sex = c("M", "F", "Both"),
  yld = c(val(3, 1, 1), val(3, 1, 2), val(3, 1, 1) + val(3, 1, 2)),
  prevalent_cases = c(val(5, 1, 1), val(5, 1, 2), val(5, 1, 1) + val(5, 1, 2)))
gbd_check$yld_per_case <- gbd_check$yld / gbd_check$prevalent_cases
print(gbd_check, digits = 4)
cat("YLD per case: GBD 2023 (anxiety disorders) =", round(gbd_check$yld_per_case[3], 4),
    "| this estimate (PTSD) =", round(q3(sim$yld_case)[2], 4), "\n")

POP23_M <- val(3, 1, 1) / val(3, 3, 1) * 1e5        # number / rate
POP23_F <- val(3, 1, 2) / val(3, 3, 2) * 1e5
cat("US population aged 25+, 2023 (GBD): M", format(round(POP23_M), big.mark = ","),
    "| F", format(round(POP23_F), big.mark = ","), "\n")
sim23 <- run_ptsd_daly(DW, SEV, PREV, POP23_M, POP23_F, seed = 18)
mc <- function(label, x) { q <- q3(x)
  data.frame(measure = label, median = q[2], low_95 = q[1], high_95 = q[3]) }
summ23 <- rbind(
  mc("PTSD cases, US adults 25+ (2023)",                     sim23$n_cases),
  mc("DALY per year, US adults 25+ (2023)",                  sim23$yld_total),
  mc("DALY per 100,000 adults 25+ per year (2023)",          sim23$rate_100k),
  data.frame(measure = "GBD 2023 check: YLD per prevalent case, anxiety disorders, 25+",
             median = gbd_check$yld_per_case[3], low_95 = NA, high_95 = NA))
summ23$axis <- AXIS; summ23$condition <- "PTSD, US adults 25+, 2023 (harmonized with axes C, D)"
print(summ23, digits = 4)
write.csv(summ23, file.path(RES_DIR, paste0(AXIS, "_PTSD_DALY_2023.csv")), row.names = FALSE)
