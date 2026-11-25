################################################################################
#      US - ITALY COMPARISON (DALY, SAME METHOD, ITALIAN/EUROPEAN INPUTS)
################################################################################
#namefile: DALY/06_italy_comparison.R
#creation dt: 14/10/2026
#
# Repeats Axes A, C and D for Italy (2023, adults aged 25+) with the same
# method as scripts 01, 03 and 04, replacing only the population inputs.
# Axis B (US military) has no Italian counterpart and is not repeated.
# No monetary values: Italy has no reference value of a DALY comparable to
# the US $100-150k (Neumann 2014); the comparison is in DALYs.
#
# Italian/European inputs:
#   A - 12-month PTSD prevalence, Italy, DSM-IV, ESEMeD 2001-02, N = 4,712:
#       0.7% (SE 0.2) [Carmassi et al. 2014, J Psychiatr Res 59:77-84,
#       doi 10.1016/j.jpsychires.2014.09.006]. Total only (no estimate by sex).
#       Weights and severity as for Axis A (GBD 2019).
#   C - individual ABSOLUTE poverty rate, Italy 2023: ages 18-34 11.8%, 35-64
#       9.4%, 65+ 6.2% [ISTAT, Le statistiche dell'Istat sulla poverta', 2023,
#       17/10/2024]. Prevalence 25-64 = mean weighted by the GBD population
#       (25-34 at 11.8%, 35-64 at 9.4%). ISTAT gives no interval, so the value
#       is fixed. RRs of Galea 2011 (US) as for Axis C.
#   D - pooled ACE RRs and prevalence for EUROPE [Bellis et al. 2019, Table 2
#       and Results]: 1 ACE 23.5% (18.7-28.5), 2+ ACEs 18.7% (14.7-23.2).
#   GBD 2023, Italy: data/source/GBD2023_*_ITA/
#
# Main limits (Chapter 3): the Italian PTSD prevalence is DSM-IV and 20 years
# old (US: DSM-5, 2012-13); ISTAT poverty is absolute and based on
# CONSUMPTION, while the RRs of Galea are estimated on INCOME thresholds in
# the US (transferability of the RRs not verified); European ACE studies are
# few (1-5 per cause).
#
# Run AFTER 01, 03 and 04 (it reads the US results):
#   source("DALY/06_italy_comparison.R", echo = TRUE)
################################################################################

source("DALY/00_daly_functions.R")
N_SIM <- 10000
one_csv <- function(dir) {
  f <- list.files(file.path(SRC_DIR, dir), "\\.csv$", full.names = TRUE)
  stopifnot(length(f) == 1); f
}
mc <- function(axis, label, x) { q <- q3(x)
  data.frame(axis = axis, measure = label, median = q[2], low_95 = q[1], high_95 = q[3]) }

# Axis A: PTSD ----
anx <- read.csv(one_csv("GBD2023_CAUSE_anxiety_ITA"), stringsAsFactors = FALSE)
stopifnot(unique(anx$year) == 2023, length(unique(anx$age_name)) == 1)
val <- function(measure_id, metric_id, sex_id)
  anx$val[anx$measure_id == measure_id & anx$metric_id == metric_id & anx$sex_id == sex_id]
POP_A <- sum(sapply(1:2, function(s) val(3, 1, s) / val(3, 3, s) * 1e5))   # M + F, aged 25+
gbd_yld_case_ita <- (val(3, 1, 1) + val(3, 1, 2)) / (val(5, 1, 1) + val(5, 1, 2))
cat("Italy, aged 25+ (GBD 2023):", format(round(POP_A), big.mark = ","),
    "| YLD per prevalent case, anxiety disorders:", round(gbd_yld_case_ita, 4), "\n")

PREV_ITA <- c(est = .007, low = .007 - 1.96 * .002, high = .007 + 1.96 * .002)   # Carmassi 2014
DW  <- read_disability_weights(c(mild = "Mild anxiety disorders", moderate = "Moderate anxiety disorders",
                                 severe = "Severe anxiety disorders"))
# a single prevalence (not by sex): the whole population in pop_m, pop_f = 0
simA <- run_ptsd_daly(DW, SEVERITY_ANXIETY, rbind(M = PREV_ITA, F = PREV_ITA), POP_A, 0, seed = 18)
resA <- rbind(mc("A", "DALY per case per year", simA$yld_case),
              mc("A", "Exposed (PTSD cases)", simA$n_cases),
              mc("A", "DALY per year", simA$yld_total),
              mc("A", "DALY per 100,000 adults 25+ per year", simA$rate_100k),
              data.frame(axis = "A", measure = "GBD 2023 check: YLD per prevalent case, anxiety disorders, Italy",
                         median = gbd_yld_case_ita, low_95 = NA, high_95 = NA))

# Axis C: poverty ----
GROUPS  <- list(`25-64` = c(25, 64), `65+` = c(65, 95))
RR_C    <- rbind(`25-64` = c(est = 1.75, low = 1.51, high = 2.04),          # Galea 2011
                 `65+`   = c(est = 1.40, low = 1.37, high = 1.43))
fileC   <- one_csv("GBD2023_CAUSE_allcauses_ITA")
gbdC    <- read_gbd_by_age(fileC, GROUPS)
stopifnot(all(gbdC$n_age_bands[gbdC$group == "25-64"] == 8), all(gbdC$n_age_bands[gbdC$group == "65+"] == 7))
popC    <- gbd_population_by_age(fileC, GROUPS)
pop2534 <- gbd_population_by_age(fileC, list(a = c(25, 34)))
PREV_C  <- c(`25-64` = (pop2534 * .118 + (popC["25-64"] - pop2534) * .094) / popC["25-64"],
             `65+`   = .062)
names(PREV_C) <- names(GROUPS)
cat("Absolute poverty rate (ISTAT 2023): 25-64", round(PREV_C[1], 4), "| 65+", PREV_C[2], "\n")
set.seed(18)
byC <- lapply(names(GROUPS), function(g) {
  get <- function(m) { r <- gbdC[gbdC$group == g & gbdC$measure == m, ]; fit_lognormal(r$est, r$low, r$high) }
  f <- paf(PREV_C[g], rlnorm_fit(N_SIM, fit_lognormal(RR_C[g, 1], RR_C[g, 2], RR_C[g, 3])))
  data.frame(deaths = f * rlnorm_fit(N_SIM, get("deaths")), yll = f * rlnorm_fit(N_SIM, get("yll")))
})
yllC <- byC[[1]]$yll + byC[[2]]$yll; dthC <- byC[[1]]$deaths + byC[[2]]$deaths
expC <- sum(PREV_C * popC)
resC <- rbind(mc("C", "Exposed (adults in absolute poverty)", rep(expC, N_SIM)),
              mc("C", "Attributable deaths per year", dthC),
              mc("C", "Attributable YLL per year", yllC),
              mc("C", "Attributable YLL per 100,000 adults 25+ per year", yllC / sum(popC) * 1e5),
              mc("C", "Attributable YLL per person in poverty per year", yllC / expC))

# Axis D: ACEs (Bellis 2019, Europe) ----
RR1 <- rbind(anxiety = c(1.44, 1.17, 1.77), depression = c(1.54, 1.51, 1.57), cancer = c(1.08, 0.89, 1.30),
             cvd = c(1.11, 0.99, 1.24), diabetes = c(1.02, 0.88, 1.19), respiratory = c(1.19, 1.04, 1.35))
RR2 <- rbind(anxiety = c(2.56, 2.19, 2.98), depression = c(2.34, 2.19, 2.50), cancer = c(1.58, 1.32, 1.91),
             cvd = c(1.57, 1.20, 2.06), diabetes = c(1.43, 1.11, 1.83), respiratory = c(1.98, 1.42, 2.76))
PREV_D <- rbind(ace1 = c(.235, .187, .285), ace2p = c(.187, .147, .232))
CAUSE_MAP <- c("Disturbi d'ansia" = "anxiety", "Disturbi depressivi" = "depression", "Neoplasie" = "cancer",
               "Malattie cardiovascolari" = "cvd", "Diabete mellito" = "diabetes",
               "Malattie respiratorie croniche" = "respiratory")
gD <- read.csv(one_csv("GBD2023_CAUSE_6causes_ITA"), stringsAsFactors = FALSE)
gD$cause <- CAUSE_MAP[gD$cause_name]; stopifnot(!any(is.na(gD$cause)))
dl  <- gD[gD$measure_id == 2 & gD$metric_id == 1, ]; dl <- dl[match(rownames(RR1), dl$cause), ]
rt  <- gD[gD$measure_id == 2 & gD$metric_id == 3, ]; rt <- rt[match(rownames(RR1), rt$cause), ]
POP_D <- median(dl$val / rt$val * 1e5)
set.seed(18)
p1 <- rbeta_fit(N_SIM, fit_beta(PREV_D[1, 1], PREV_D[1, 2], PREV_D[1, 3]))
p2 <- rbeta_fit(N_SIM, fit_beta(PREV_D[2, 1], PREV_D[2, 2], PREV_D[2, 3]))
dalyD <- Reduce(`+`, lapply(seq_len(nrow(RR1)), function(i) {
  r1 <- rlnorm_fit(N_SIM, fit_lognormal(RR1[i, 1], RR1[i, 2], RR1[i, 3]))
  r2 <- rlnorm_fit(N_SIM, fit_lognormal(RR2[i, 1], RR2[i, 2], RR2[i, 3]))
  x  <- p1 * (r1 - 1) + p2 * (r2 - 1)
  x / (1 + x) * rlnorm_fit(N_SIM, fit_lognormal(dl$val[i], dl$lower[i], dl$upper[i]))
}))
rk <- read.csv(one_csv("GBD2023_RISK_child_maltreatment_ITA"), stringsAsFactors = FALSE)
rk <- rk[rk$metric_id == 1 & rk$measure_id == 2 & rk$rei_name == "Maltrattamento di minori", ]
resD <- rbind(mc("D", "Exposed (adults with >=1 ACE)", (p1 + p2) * POP_D),
              mc("D", "Attributable DALY per year, 6 causes", dalyD),
              mc("D", "Attributable DALY per 100,000 adults 25+ per year", dalyD / POP_D * 1e5),
              mc("D", "Attributable DALY per person with >=1 ACE per year", dalyD / ((p1 + p2) * POP_D)),
              data.frame(axis = "D", measure = "Check: GBD 2023 DALY attributable to child maltreatment, Italy",
                         median = rk$val, low_95 = rk$lower, high_95 = rk$upper))

ita <- rbind(resA, resC, resD); ita$country <- "Italy"; rownames(ita) <- NULL
print(ita, digits = 4)
write.csv(ita, file.path(RES_DIR, "ITA_DALY.csv"), row.names = FALSE)

# Comparison with the US (same measures, 2023, adults aged 25+) ----
rd <- function(f) read.csv(file.path(RES_DIR, f), stringsAsFactors = FALSE)
pick <- function(d, m) { r <- d[d$measure == m, ]; stopifnot(nrow(r) == 1); unlist(r[, c("median", "low_95", "high_95")]) }
A <- rd("A_PTSD_DALY.csv"); A23 <- rd("A_PTSD_DALY_2023.csv"); C <- rd("C_poverty_DALY.csv"); D <- rd("D_ACE_DALY.csv")
pi <- function(m) pick(ita, m)
cmp <- rbind(
  data.frame(axis = "A", metric = "per person exposed", country = c("United States", "Italy"),
             rbind(pick(A, "DALY per case per year"), pi("DALY per case per year"))),
  data.frame(axis = "A", metric = "per 100,000 adults", country = c("United States", "Italy"),
             rbind(pick(A23, "DALY per 100,000 adults 25+ per year (2023)"), pi("DALY per 100,000 adults 25+ per year"))),
  data.frame(axis = "C", metric = "per person exposed", country = c("United States", "Italy"),
             rbind(pick(C, "Attributable YLL per person in poverty per year"), pi("Attributable YLL per person in poverty per year"))),
  data.frame(axis = "C", metric = "per 100,000 adults", country = c("United States", "Italy"),
             rbind(pick(C, "Attributable YLL per 100,000 adults 25+ per year"), pi("Attributable YLL per 100,000 adults 25+ per year"))),
  data.frame(axis = "D", metric = "per person exposed", country = c("United States", "Italy"),
             rbind(pick(D, "Attributable DALY per person with >=1 ACE per year"), pi("Attributable DALY per person with >=1 ACE per year"))),
  data.frame(axis = "D", metric = "per 100,000 adults", country = c("United States", "Italy"),
             rbind(pick(D, "Attributable DALY per 100,000 adults 25+ per year"), pi("Attributable DALY per 100,000 adults 25+ per year"))))
rownames(cmp) <- NULL
print(cmp, digits = 3)
write.csv(cmp, file.path(RES_DIR, "ITA_vs_USA.csv"), row.names = FALSE)

lab <- c(A = "A - PTSD", C = "C - Poverty (YLL)", D = "D - Childhood adversity")
cmp$label  <- factor(lab[cmp$axis], levels = rev(lab))
cmp$metric <- factor(cmp$metric, levels = c("per person exposed", "per 100,000 adults"),
                     labels = c("Per person exposed per year", "Per 100,000 adults per year"))
cmp$country <- factor(cmp$country, levels = c("United States", "Italy"))
g <- ggplot(cmp, aes(median, label, colour = country, shape = country)) +
  geom_errorbar(aes(xmin = low_95, xmax = high_95), orientation = "y", width = 0.25,
                linewidth = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_colour_manual(values = c("United States" = "#2a78d6", "Italy" = "#eb6834"), name = NULL) +
  scale_shape_manual(values = c("United States" = 16, "Italy" = 17), name = NULL) +
  expand_limits(x = 0) +
  labs(x = "Years of healthy life lost (median and 95% interval)", y = NULL,
       title = "Annual health burden by axis: United States vs Italy, adults 25+, 2023",
       caption = "Same method; country-specific population inputs. Axis B (US military) has no Italian equivalent.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"), plot.caption = element_text(size = 8, colour = "grey30"))
ggsave(file.path(FIG_DIR, "X_USA_vs_ITA.png"), g, width = 9, height = 4.2, dpi = 300)
ggsave(file.path(FIG_DIR, "X_USA_vs_ITA.pdf"), g, width = 9, height = 4.2)
