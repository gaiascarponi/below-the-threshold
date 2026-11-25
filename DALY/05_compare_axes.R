################################################################################
#    COMPARISON ACROSS AXES - BURDEN OF DISEASE (SIDE BY SIDE, NOT SUMMED)
################################################################################
#namefile: DALY/05_compare_axes.R
#creation dt: 13/10/2026
#
# Reads the results of the four axis scripts and produces a comparison table
# and figure. The estimates are NOT summed (Chapter 3, "Interpretation and
# non-additivity"): A and B are the burden of a disease (YLD only), C the
# burden attributable to an exposure (YLL only), D the burden attributable to
# an exposure (YLL + YLD). The DALY component is shown in the table and in the
# figure (color and marker shape).
#
# Axis A: the harmonized 2023 estimate is used (adults aged 25+, GBD 2023
# population; last section of 01_axisA), on the same basis as C and D. The
# DALY per case is the same.
# Axis B: US military population of Davis 2022 (not available in GBD).
#
# Two measures, in two separate panels (different scales, never a double axis):
#   - per person exposed per year (case of PTSD / person in poverty / person
#     with 1+ ACE), also expressed in days of healthy life
#   - per 100,000 adults per year
#
# Run AFTER scripts 01-04:
#   source("DALY/05_compare_axes.R", echo = TRUE)
################################################################################

source("DALY/00_daly_functions.R")

rd <- function(f) read.csv(file.path(RES_DIR, f), stringsAsFactors = FALSE)
pick <- function(d, m) { r <- d[d$measure == m, ]; stopifnot(nrow(r) == 1)
  c(median = r$median, low = r$low_95, high = r$high_95) }

A <- rd("A_PTSD_DALY.csv"); A23 <- rd("A_PTSD_DALY_2023.csv"); B <- rd("B_PTSD_DALY.csv")
C <- rd("C_poverty_DALY.csv"); D <- rd("D_ACE_DALY.csv")

axes <- data.frame(
  axis      = c("A", "B", "C", "D"),
  label     = c("A - PTSD, adults", "B - PTSD, military",
                "C - Poverty", "D - Childhood adversity (ACEs)"),
  exposure  = c("case of PTSD", "case of PTSD", "person in poverty", "person with >=1 ACE"),
  type      = c("disease", "disease", "attributable", "attributable"),
  component = c("YLD only", "YLD only", "YLL only", "YLL + YLD"))

per_exp <- rbind(pick(A, "DALY per case per year"),
                 pick(B, "DALY per case per year"),
                 pick(C, "Attributable YLL per person in poverty per year"),
                 pick(D, "Attributable DALY per person with >=1 ACE per year"))
per_100k <- rbind(pick(A23, "DALY per 100,000 adults 25+ per year (2023)"),
                  pick(B, "DALY per 100,000 adults per year"),
                  pick(C, "Attributable YLL per 100,000 adults 25+ per year"),
                  pick(D, "Attributable DALY per 100,000 adults 25+ per year"))
total <- rbind(pick(A23, "DALY per year, US adults 25+ (2023)"), pick(B, "DALY per year"),
               pick(C, "Attributable YLL per year, 25+"),
               pick(D, "Attributable DALY per year, 6 causes"))
usd <- rbind(pick(A, "Monetized DALY per case, USD 2024"),
             pick(B, "Monetized DALY per case, USD 2024"),
             pick(C, "Monetized YLL per person in poverty, USD 2024"),
             pick(D, "Monetized DALY per person with >=1 ACE, USD 2024"))

tab <- cbind(axes,
             per_exposed_median = per_exp[, "median"], per_exposed_low = per_exp[, "low"],
             per_exposed_high = per_exp[, "high"],
             days_per_exposed = per_exp[, "median"] * 365.25,
             per_100k_median = per_100k[, "median"], per_100k_low = per_100k[, "low"],
             per_100k_high = per_100k[, "high"],
             total_median = total[, "median"], total_low = total[, "low"], total_high = total[, "high"],
             usd_per_exposed_median = usd[, "median"])
print(tab, digits = 3)
write.csv(tab, file.path(RES_DIR, "X_axes_comparison.csv"), row.names = FALSE)

# Figure: point (median) and 95% interval, two independent panels ----
long <- rbind(
  data.frame(tab[, c("label", "component")], panel = "Per person exposed per year",
             median = tab$per_exposed_median, low = tab$per_exposed_low, high = tab$per_exposed_high),
  data.frame(tab[, c("label", "component")], panel = "Per 100,000 adults per year",
             median = tab$per_100k_median, low = tab$per_100k_low, high = tab$per_100k_high))
long$panel <- factor(long$panel, levels = c("Per person exposed per year", "Per 100,000 adults per year"))
long$label <- factor(long$label, levels = rev(axes$label))
long$component <- factor(long$component, levels = c("YLD only", "YLL only", "YLL + YLD"))

# validated categorical palette; color is always paired with the marker shape
PAL <- c("YLD only" = "#2a78d6", "YLL only" = "#eb6834", "YLL + YLD" = "#1baf7a")
SHP <- c("YLD only" = 16, "YLL only" = 17, "YLL + YLD" = 15)

g <- ggplot(long, aes(x = median, y = label, colour = component, shape = component)) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", height = 0.18, linewidth = 0.6) +
  geom_point(size = 3) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_colour_manual(values = PAL, name = "DALY component") +
  scale_shape_manual(values = SHP, name = "DALY component") +
  expand_limits(x = 0) +
  labs(x = "Years of healthy life lost (median and 95% interval)", y = NULL,
       title = "Annual health burden by axis, United States",
       caption = "Estimates are compared side by side and are not additive. A, B: disease burden; C, D: burden attributable to an exposure.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
        plot.caption = element_text(size = 8, colour = "grey30"))
ggsave(file.path(FIG_DIR, "X_axes_comparison.png"), g, width = 9, height = 4.2, dpi = 300)
ggsave(file.path(FIG_DIR, "X_axes_comparison.pdf"), g, width = 9, height = 4.2)
