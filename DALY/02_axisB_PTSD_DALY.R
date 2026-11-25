################################################################################
#                 AXIS B - BURDEN OF PTSD (DALY) - US MILITARY
################################################################################
#namefile: DALY/02_axisB_PTSD_DALY.R
#creation dt: 13/10/2026
#
# Social cost model, separate from the transcriptomic analysis: it does not use
# the expression data. It estimates the burden of PTSD (the clinical outcome of
# Axis B, GSE64813, US Marines) in the US military population (active duty and
# veterans, as in Davis 2022), in DALYs, with the GBD method.
# Same model as Axis A (01_axisA_PTSD_DALY.R): only prevalences, populations
# and costs change. The DALY per case is therefore equal to Axis A by
# construction: GBD does not distinguish the severity of combat PTSD (a limit).
#
#   DALY = YLL + YLD;  YLL = 0 (GBD assigns no direct mortality to anxiety
#   disorders);  YLD per case = sum_s p_s * DW_s
#   PTSD is not a separate GBD cause: it is part of "anxiety disorders"
#   (ICD-10 F43.1, DSM-IV-TR 309.81; Yang et al. 2021).
#
# Inputs:
#   - GBD 2019 disability weights (IHME 2020), file in data/source
#   - GBD 2019 severity distribution of anxiety disorders (Yang et al. 2021,
#     Epidemiol Psychiatr Sci 30:e36, supplementary) -> 00_daly_functions.R
#   - 12-month PTSD prevalence, VETERANS (DSM-5, NESARC-III 2012-13): M 6.7%
#     (5.7-7.8), F 11.7% (7.1-16.4) [Lehavot et al. 2018, Am J Prev Med
#     54:e1-e9]; Davis 2022 applies it to the whole military population
#     (active duty and veterans). The interval for women is wide (few women).
#   - US military populations 2018 and costs per case (USD 2018) from Davis
#     et al. 2022, J Clin Psychiatry 83(3):21m14116: 1,368,925 M and 293,256 F
#     with PTSD; total societal cost $25,684; direct health care costs $6,071
#     (Table 1, military health care and VA)
#   - CPI-U CUUR0000SA0 (BLS), file in data/source
#
# Uncertainty: as for Axis A (beta by quantile matching, uniform for the value
# of a DALY, 10,000 draws, seed 18; checks with seed 2026 and Dirichlet).
#
# Outputs:
#   results/B_PTSD_DALY.csv               main estimates
#   results/B_distribution_fit.csv        published vs reproduced (beta / normal)
#   results/B_sensitivity.csv             stability and sensitivity
#   figures/B_PTSD_DALY_vs_cost.png/.pdf
#
# Run from the project root:
#   source("DALY/02_axisB_PTSD_DALY.R", echo = TRUE)
################################################################################

source("DALY/00_daly_functions.R")

AXIS      <- "B"
CONDITION <- "PTSD (US military: active duty + veterans)"
POP_M     <- 1368925 / 0.067          # US military population 2018, men (Davis 2022)
POP_F     <- 293256 / 0.117           # US military population 2018, women (Davis 2022)
PREV      <- rbind(M = c(est = .067, low = .057, high = .078),     # Lehavot 2018, veterans
                   F = c(est = .117, low = .071, high = .164))
COST_TOT  <- 25684                    # USD 2018 per case, total societal cost
COST_HC   <- 6071                     # USD 2018 per case, direct health care costs
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
plot_ptsd_costs(summ, "Axis B - PTSD (US military): annual burden per case",
                paste0(AXIS, "_PTSD_DALY_vs_cost"))
