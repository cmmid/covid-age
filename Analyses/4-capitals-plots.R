library(Rcpp) 
library(RcppGSL)
library(HDInterval)
library(ggplot2)
library(data.table)
library(socialmixr)
library(lubridate)
library(stringr)
library(Hmisc)
library(extraDistr)
library(countrycode)
library(nloptr)
library(qs)
library(rlang)
library(readxl)
library(cowplot)

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
cm_verbose = F
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_builddir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# Return the mean and CI from a set of samples x.
mean_ci = function(x, ci = 0.95)
{
    h = quantile(x, c(0.5 - ci/2, 0.5 + ci/2));
    m = mean(x);
    return (list(mean = m, lower = h[[1]], upper = h[[2]]))
}

# Load and process data
cap = cm_load(path("4-capitals-fix-u.qs"))
gbd = fread("~/Dropbox/nCoV/gbd.regions.txt");
gbd$iso2c = countrycode(gbd$Country, origin = "country.name", dest = "iso2c")

gbd[SuperRegion %like% "Latin", SuperRegionA := "LAC"]
gbd[SuperRegion %like% "High", SuperRegionA := "HI"]
gbd[SuperRegion %like% "Central", SuperRegionA := "CE, EE, CA"]
gbd[SuperRegion %like% "Sub", SuperRegionA := "SSA"]
gbd[SuperRegion %like% "Southeast", SuperRegionA := "SEA, EA, O"]
gbd[SuperRegion %like% "Middle", SuperRegionA := "NA, ME"]
gbd[SuperRegion %like% "South Asia", SuperRegionA := "SA"]

epi = cap$epi[virus == "cov" & compartment == "cases"]
dist = cap$dist
age_info = cap$age_info
mean_age = cap$mean_age

#Plotting
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# B. Wiggly plot
res2 = cap$epi[compartment == "cases" & virus == "cov" & schools == "open" & fIa == 0.5]
pd = res2[, .(peak_day = which.max(cases) - 1), by = .(run, location)]
res2 = merge(res2, pd, by = c("run", "location"))
ma2 = mean_age[what == "cases" & virus == "cov" & fIa == 0.5]
ma2 = merge(ma2, pd, by = c("run", "location"))

res2 = merge(res2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
res2 = merge(res2, gbd, by = "iso2c")

ma2 = merge(ma2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
ma2 = merge(ma2, gbd, by = "iso2c")

res2h = res2[, cm_mean_hdi(cases/popK), by = .(location, t_minus_peak = t - peak_day, SuperRegionA)]
ma2h = ma2[, cm_mean_hdi(mean_age), by = .(location, t_minus_peak = t - peak_day, SuperRegionA)]

f4b = ggplot() +
    geom_ribbon(data = ma2h[t_minus_peak < 200], aes(t_minus_peak, ymin = lower, ymax = upper, group = location, fill = SuperRegionA), alpha = 0.1) +
    geom_line(data = ma2h[t_minus_peak < 200], aes(t_minus_peak, mean, group = location, colour = SuperRegionA), alpha = 0.5, size = 0.25) +
    geom_ribbon(data = res2h[t_minus_peak < 200], aes(t_minus_peak, ymin = lower, ymax = upper, group = location, fill = SuperRegionA), alpha = 0.1) +
    geom_line(data = res2h[t_minus_peak < 200], aes(t_minus_peak, mean, group = location, colour = SuperRegionA), alpha = 0.5, size = 0.25) +
    theme(legend.position = "none") +
    labs(x = "Time relative to peak (days)", y = "Mean age of cases +\nClinical cases per 1000 population")

ggsave(path("../Submission/Fig4b.pdf"), f4b, width = 10, height = 12, units = "cm", useDingbats = F)


# C. Distribution of cases by country / region
dist = merge(dist, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
dist = merge(dist, gbd, by = "iso2c")
dist[, age_mid := as.numeric(str_sub(age, 1, 2)) + 5]

dist[period == "early", period := "Early"]
dist[period == "middle", period := "Middle"]
dist[period == "late", period := "Late"]

dist = dist[, cm_mean_hdi(cases/popK), by = .(virus, fIa, period, age_mid, country, SuperRegionA)]

f4c = ggplot(dist[virus == "cov" & fIa == 0.5 & period != "Middle"]) + 
    geom_line(aes(x = age_mid, y = mean, group = country, colour = SuperRegionA), alpha = 0.2, size = 0.25) +
    facet_grid(SuperRegionA~period) +
    theme(legend.position = "none", strip.background = element_blank()) +
    labs(x = "Age", y = "Clinical cases per 1000 population") +
    xlim(0, 75)

ggsave(path("../Submission/Fig4c.pdf"), f4c, width = 7, height = 12, units = "cm", useDingbats = F)

# f4d = ggplot(dist2[virus == "cov" & fIa == 0.5]) + 
#     geom_line(aes(x = age_mid, y = incidence, colour = period, group = period), alpha = 0.5) +
#     scale_y_log10() +
#     facet_wrap(~SuperRegion, scales = "free", ncol = 1) +
#     theme(legend.position = "bottom")
# 
# ggsave(path("../Submission/Fig4d.pdf"), f4c, width = 5, height = 12, units = "cm", useDingbats = F)

# A. Scatter plot
epi = cap$epi[virus == "cov" & compartment == "cases"]
epi[, week := t %/% 7]
epi = epi[, .(incidence = sum(cases)), by = .(run, virus, location, compartment, fIa, week, popK)]
epi2 = epi[, .(peak_incidence = max(incidence), total_cases = sum(incidence)), by = .(run, virus, location, compartment, fIa, popK)]

epi2 = merge(epi2, age_info, by.x = "location", by.y = "name")
epi2 = merge(epi2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
epi2 = merge(epi2, gbd, by = "iso2c")

epi2p = epi2[, cm_mean_hdi(peak_incidence / popK), by = .(median_age, virus, location, compartment, fIa, popK, SuperRegion)]
epi2t = epi2[,    cm_mean_hdi(total_cases / popK), by = .(median_age, virus, location, compartment, fIa, popK, SuperRegion)]

# Version with 3 panels
f4Sa = ggplot() + 
    geom_smooth(data = epi2p, aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange(data = epi2p, aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point(data = epi2p, aes(median_age, mean, shape = "peak", colour = SuperRegion), alpha = 0.75, size = 2) +
    geom_smooth(data = epi2t, aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange(data = epi2t, aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point(data = epi2t, aes(median_age, mean, shape = "total", colour = SuperRegion), alpha = 0.75, size = 2) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17)) +
    labs(x = "Median age", y = "Attack rate per 1000 population", colour = "Region", shape = "Incidence") +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    facet_wrap(~fIa) +
    theme(legend.position = "bottom") 
f4Sa
ggsave(path("../Submission/FigS_countries_fIa_clinical_infections.pdf"), f4Sa, width = 18, height = 18, units = "cm", useDingbats = F)

f4a = ggplot() + 
    geom_smooth     (data = epi2p[fIa == 0.5], aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange  (data = epi2p[fIa == 0.5], aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point      (data = epi2p[fIa == 0.5], aes(median_age, mean, shape = "peak", colour = SuperRegion), alpha = 0.75, size = 2) +
    geom_smooth     (data = epi2t[fIa == 0.5], aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange  (data = epi2t[fIa == 0.5], aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point      (data = epi2t[fIa == 0.5], aes(median_age, mean, shape = "total", colour = SuperRegion), alpha = 0.75, size = 2) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17), breaks = c("total", "peak"), labels = c("Total", "Peak")) +
    labs(x = "Median age", y = "Clinical attack rate per 1000 population", colour = NULL, shape = NULL) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) + txt_theme +
    theme(legend.position = c(0.03, 0.888), legend.key.height = unit(0.15, "cm"), legend.box = "horizontal") 

ggsave(path("../Submission/Fig4a.pdf"), f4a, width = 10, height = 10, units = "cm", useDingbats = F)

# Subclinical infections version
epi = cap$epi[virus == "cov" & compartment == "subclinical"]
epi[, week := t %/% 7]
epi = epi[, .(incidence = sum(cases)), by = .(run, virus, location, compartment, fIa, week, popK)]
epi2 = epi[, .(peak_incidence = max(incidence), total_cases = sum(incidence)), by = .(run, virus, location, compartment, fIa, popK)]

epi2 = merge(epi2, age_info, by.x = "location", by.y = "name")
epi2 = merge(epi2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
epi2 = merge(epi2, gbd, by = "iso2c")

epi2p = epi2[, cm_mean_hdi(peak_incidence / popK), by = .(median_age, virus, location, compartment, fIa, popK, SuperRegion)]
epi2t = epi2[,    cm_mean_hdi(total_cases / popK), by = .(median_age, virus, location, compartment, fIa, popK, SuperRegion)]

f4Sai = ggplot() + 
    geom_smooth(data = epi2p, aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange(data = epi2p, aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point(data = epi2p, aes(median_age, mean, shape = "peak", colour = SuperRegion), alpha = 0.75, size = 2) +
    geom_smooth(data = epi2t, aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange(data = epi2t, aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point(data = epi2t, aes(median_age, mean, shape = "total", colour = SuperRegion), alpha = 0.75, size = 2) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17)) +
    labs(x = "Median age", y = "Subclinical attack rate per 1000 population", colour = "Region", shape = "Incidence") +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    facet_wrap(~fIa) +
    theme(legend.position = "bottom") 

ggsave(path("../Submission/FigS_countries_fIa_subclinical_infections.pdf"), f4Sai, width = 18, height = 18, units = "cm", useDingbats = F)

f4ai = ggplot() + 
    geom_smooth     (data = epi2p[fIa == 0.5], aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange  (data = epi2p[fIa == 0.5], aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point      (data = epi2p[fIa == 0.5], aes(median_age, mean, shape = "peak", colour = SuperRegion), alpha = 0.75, size = 1) +
    geom_smooth     (data = epi2t[fIa == 0.5], aes(median_age, mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange  (data = epi2t[fIa == 0.5], aes(median_age, ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point      (data = epi2t[fIa == 0.5], aes(median_age, mean, shape = "total", colour = SuperRegion), alpha = 0.75, size = 1) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17), breaks = c("total", "peak"), labels = c("Total", "Peak")) +
    labs(x = "Median age", y = "Subclinical attack rate\nper 1000 population", colour = NULL, shape = NULL) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) + txt_theme +
    theme(legend.position = "none") 

ggsave(path("../Submission/FigS_4a_subclinical_infections.pdf"), f4ai, width = 10, height = 10, units = "cm", useDingbats = F)


# R0
R0s = cap$R0
R0s = merge(R0s, age_info, by.x = "location", by.y = "name")
R0s = merge(R0s, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
R0s = merge(R0s, gbd, by = "iso2c")

R0s = R0s[, cm_mean_hdi(R0), by = .(median_age, SuperRegion, fIa)]

# by fIa
fR0f = ggplot(R0s, aes(x = median_age)) + 
    geom_smooth(aes(y = mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.25, size = 0.25) +
    geom_point(aes(y = mean, colour = SuperRegion), alpha = 0.75) +
    labs(x = "Median age", y = expression(R[0])) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    facet_wrap(~paste("f_Ia =", fIa)) +
    theme(legend.position = "none", strip.background = element_blank()) 

ggsave(path("../Submission/FigS_R0s.pdf"), fR0f, width = 12, height = 5, units = "cm", useDingbats = F)

# Just 0.5
fR0 = ggplot(R0s[fIa == 0.5], aes(x = median_age)) + 
    geom_smooth(aes(y = mean), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper, colour = SuperRegion), alpha = 0.5, size = 0.5) +
    geom_point(aes(y = mean, colour = SuperRegion), alpha = 0.75, size = .75) +
    labs(x = "Median age", y = expression(R[0])) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    theme(legend.position = "none", strip.background = element_blank()) 

ggsave(path("../Submission/Fig4c.pdf"), fR0, width = 6, height = 6, units = "cm", useDingbats = F)

f4 = plot_grid(
    f4a,
    plot_grid(f4ai, fR0, f4b, ncol = 1, labels = c("b", "c", "d"), label_size = 9, label_x = -0.05, align = "hv", axis = "b"),
    f4c, 
    nrow = 1, rel_widths = c(5, 3, 4), labels = c("a", "", "e"), label_size = 9)

ggsave(path("../Submission/Fig4.pdf"), f4, width = 22, height = 10, units = "cm", useDingbats = F)









# VERSION FOR LOW INCOME COUNTRIES
if (0) {

library(Rcpp) 
library(RcppGSL)
library(HDInterval)
library(ggplot2)
library(data.table)
library(socialmixr)
library(lubridate)
library(stringr)
library(Hmisc)
library(extraDistr)
library(countrycode)
library(nloptr)
library(qs)
library(rlang)
library(readxl)
library(cowplot)
    
# worldbank income
worldbank = fread("~/Dropbox/nCoV/WorldBankIncome.csv")
worldbank[, iso2c := countrycode(Code, "wb", "iso2c")]
worldbank[, income := "Higher"]
worldbank[`Income group` %in% c("Lower middle income", "Low income"), income := "Lower"]

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
cm_verbose = F
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_builddir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# Return the mean and CI from a set of samples x.
mean_ci = function(x, ci = 0.95)
{
    h = quantile(x, c(0.5 - ci/2, 0.5 + ci/2));
    m = mean(x);
    return (list(mean = m, lower = h[[1]], upper = h[[2]]))
}

# Load and process data
cap = cm_load(path("4-capitals-fix-u-lowincome.qs"))
gbd = fread("~/Dropbox/nCoV/gbd.regions.txt");
gbd$iso2c = countrycode(gbd$Country, origin = "country.name", dest = "iso2c")

gbd[SuperRegion %like% "Latin", SuperRegionA := "LAC"]
gbd[SuperRegion %like% "High", SuperRegionA := "HI"]
gbd[SuperRegion %like% "Central", SuperRegionA := "CE, EE, CA"]
gbd[SuperRegion %like% "Sub", SuperRegionA := "SSA"]
gbd[SuperRegion %like% "Southeast", SuperRegionA := "SEA, EA, O"]
gbd[SuperRegion %like% "Middle", SuperRegionA := "NA, ME"]
gbd[SuperRegion %like% "South Asia", SuperRegionA := "SA"]

epi = cap$epi[virus == "cov" & compartment == "cases"]
dist = cap$dist
age_info = cap$age_info
mean_age = cap$mean_age

#Plotting
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# B. Wiggly plot
res2 = cap$epi[compartment == "cases" & virus == "cov" & schools == "open" & fIa == 0.5]
pd = res2[, .(peak_day = which.max(cases) - 1), by = location]
res2 = merge(res2, pd, by = "location")
ma2 = mean_age[what == "cases" & virus == "cov" & fIa == 0.5]
ma2 = merge(ma2, pd, by = "location")

res2 = merge(res2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
res2 = merge(res2, gbd, by = "iso2c")

ma2 = merge(ma2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
ma2 = merge(ma2, gbd, by = "iso2c")

f4b = ggplot() +
    geom_line(data = ma2[t - peak_day < 200], aes(x = t - peak_day, y = mean_age, colour = SuperRegionA, group = location), alpha = 0.5, size = 0.25) +
    geom_line(data = res2[t - peak_day < 200], aes(t - peak_day, cases/popK, group = location, colour = SuperRegionA), alpha = 0.5, size = 0.25) +
    theme(legend.position = "none") +
    labs(x = "Time relative to peak (days)", y = "Mean age of cases +\nClinical cases per 1000 population")

ggsave(path("../Submission/Supp Figs/LowIncomeFig4b.pdf"), f4a, width = 10, height = 12, units = "cm", useDingbats = F)


# C. Distribution of cases by country / region
dist = merge(dist, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
dist = merge(dist, gbd, by = "iso2c")
dist[, age_mid := as.numeric(str_sub(age, 1, 2)) + 5]

dist[period == "early", period := "Early"]
dist[period == "middle", period := "Middle"]
dist[period == "late", period := "Late"]

dist = merge(dist, worldbank, by = "iso2c", all.x = T)

f4c = ggplot(dist[virus == "cov" & fIa == 0.5 & period != "Middle"]) + 
    geom_line(aes(x = age_mid, y = cases/popK, group = country, colour = SuperRegionA), alpha = 0.2, size = 0.25) +
    facet_grid(SuperRegionA~period) +
    theme(legend.position = "none", strip.background = element_blank()) +
    labs(x = "Age", y = "Clinical cases per 1000 population") +
    xlim(0, 75)

ggsave(path("../Submission/Supp Figs/LowIncomeFig4c.pdf"), f4c, width = 7, height = 12, units = "cm", useDingbats = F)

# f4d = ggplot(dist2[virus == "cov" & fIa == 0.5]) + 
#     geom_line(aes(x = age_mid, y = incidence, colour = period, group = period), alpha = 0.5) +
#     scale_y_log10() +
#     facet_wrap(~SuperRegion, scales = "free", ncol = 1) +
#     theme(legend.position = "bottom")
# 
# ggsave(path("../Submission/Supp Figs/LowIncomeFig4d.pdf"), f4c, width = 5, height = 12, units = "cm", useDingbats = F)

# A. Scatter plot
epi = cap$epi[virus == "cov" & compartment == "cases"]
epi[, week := t %/% 7]
epi = epi[, .(incidence = sum(cases)), by = .(virus, location, compartment, fIa, week, popK)]
epi2 = epi[, .(peak_incidence = max(incidence), total_cases = sum(incidence)), by = .(virus, location, compartment, fIa, popK)]

epi2 = merge(epi2, age_info, by.x = "location", by.y = "name")
epi2 = merge(epi2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
epi2 = merge(epi2, gbd, by = "iso2c")
epi2 = merge(epi2, worldbank, by = "iso2c", all.x = T)

# Version with 3 panels
f4Sa = ggplot(epi2, aes(x = median_age)) + 
    geom_smooth(aes(y = peak_incidence / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = peak_incidence / popK, shape = "peak", colour = SuperRegion), alpha = 0.75) +
    geom_smooth(aes(y = total_cases / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = total_cases / popK, shape = "total", colour = SuperRegion), alpha = 0.75) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17)) +
    labs(x = "Median age", y = "Attack rate per 1000 population", colour = "Region", shape = "Incidence") +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    facet_wrap(~fIa) +
    theme(legend.position = "bottom") 

ggsave(path("../Submission/Supp Figs/LowIncomeFigS_countries_fIa_clinical_infections.pdf"), f4Sa, width = 18, height = 10, units = "cm", useDingbats = F)

f4a = ggplot(epi2[fIa == 0.5], aes(x = median_age)) + 
    geom_smooth(aes(y = peak_incidence / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = peak_incidence / popK, shape = "peak", colour = SuperRegion), alpha = 0.75) +
    geom_smooth(aes(y = total_cases / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = total_cases / popK, shape = "total", colour = SuperRegion), alpha = 0.75) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17), breaks = c("total", "peak"), labels = c("Total", "Peak")) +
    labs(x = "Median age", y = "Clinical attack rate per 1000 population", colour = NULL, shape = NULL) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, 600) + txt_theme +
    theme(legend.position = c(0.03, 0.888), legend.key.height = unit(0.15, "cm"), legend.box = "horizontal") 

ggsave(path("../Submission/Supp Figs/LowIncomeFig4a.pdf"), f4a, width = 10, height = 10, units = "cm", useDingbats = F)

# Subclinical infections version
epi = cap$epi[virus == "cov" & compartment == "subclinical"]
epi[, week := t %/% 7]
epi = epi[, .(incidence = sum(cases)), by = .(virus, location, fIa, week, popK)]
epi2 = epi[, .(peak_incidence = max(incidence), total_cases = sum(incidence)), by = .(virus, location, fIa, popK)]

epi2 = merge(epi2, age_info, by.x = "location", by.y = "name")
epi2 = merge(epi2, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
epi2 = merge(epi2, gbd, by = "iso2c")
epi2 = merge(epi2, worldbank, by = "iso2c", all.x = T)

f4Sai = ggplot(epi2, aes(x = median_age)) + 
    geom_smooth(aes(y = peak_incidence / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = peak_incidence / popK, shape = "peak", colour = SuperRegion), alpha = 0.75) +
    geom_smooth(aes(y = total_cases / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = total_cases / popK, shape = "total", colour = SuperRegion), alpha = 0.75) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17)) +
    labs(x = "Median age", y = "Subclinical attack rate per 1000 population", colour = "Region", shape = "Incidence") +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    facet_wrap(~fIa) +
    theme(legend.position = "bottom") 

ggsave(path("../Submission/Supp Figs/LowIncomeFigS_countries_fIa_subclinical_infections.pdf"), f4Sai, width = 18, height = 10, units = "cm", useDingbats = F)

f4ai = ggplot(epi2[fIa == 0.5], aes(x = median_age)) + 
    geom_smooth(aes(y = peak_incidence / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = peak_incidence / popK, shape = "peak", colour = SuperRegion), alpha = 0.75) +
    geom_smooth(aes(y = total_cases / popK, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = .5) +
    geom_point(aes(y = total_cases / popK, shape = "total", colour = SuperRegion), alpha = 0.75) +
    scale_shape_manual(values = c("total" = 16, "peak" = 17), breaks = c("total", "peak"), labels = c("Total", "Peak")) +
    labs(x = "Median age", y = "Subclinical attack rate\nper 1000 population", colour = NULL, shape = NULL) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) + txt_theme +
    theme(legend.position = "none") 

ggsave(path("../Submission/Supp Figs/LowIncomeFigS_4a_subclinical_infections.pdf"), f4ai, width = 10, height = 10, units = "cm", useDingbats = F)


# R0
R0s = cap$R0
R0s = merge(R0s, age_info, by.x = "location", by.y = "name")
R0s = merge(R0s, cap$wc[, .(name, country = country.etc, iso2c = cc)], by.x = "location", by.y = "name")
R0s = merge(R0s, gbd, by = "iso2c")
R0s = merge(R0s, worldbank, by = "iso2c", all.x = T)

# by fIa
fR0f = ggplot(R0s, aes(x = median_age)) + 
    geom_smooth(aes(y = R0, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_point(aes(y = R0, colour = SuperRegion), alpha = 0.75) +
    labs(x = "Median age", y = expression(R[0])) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    facet_wrap(~paste("f_Ia =", fIa)) +
    theme(legend.position = "none", strip.background = element_blank()) 

ggsave(path("../Submission/Supp Figs/LowIncomeFigS_R0s.pdf"), fR0f, width = 12, height = 5, units = "cm", useDingbats = F)

# Just 0.5
fR0 = ggplot(R0s[fIa == 0.5], aes(x = median_age)) + 
    geom_smooth(aes(y = R0, group = income, linetype = income), method = "lm", colour = "black", alpha = 0.5, size = 0.5) +
    geom_point(aes(y = R0, colour = SuperRegion), alpha = 0.75) +
    labs(x = "Median age", y = expression(R[0])) +
    guides(shape = guide_legend(nrow = 7, byrow = T), colour = guide_legend(nrow = 7, byrow = T)) +
    ylim(0, NA) +
    theme(legend.position = "none", strip.background = element_blank()) 

ggsave(path("../Submission/Supp Figs/LowIncomeFig4c.pdf"), fR0, width = 6, height = 6, units = "cm", useDingbats = F)

f4 = plot_grid(
    f4a,
    plot_grid(f4ai, fR0, f4b, ncol = 1, labels = c("b", "c", "d"), label_size = 9, label_x = -0.05, align = "hv", axis = "b"),
    f4c, 
    nrow = 1, rel_widths = c(5, 3, 4), labels = c("a", "", "e"), label_size = 9)

ggsave(path("../Submission/Supp Figs/LowIncomeFig4.pdf"), f4, width = 22, height = 10, units = "cm", useDingbats = F)
}