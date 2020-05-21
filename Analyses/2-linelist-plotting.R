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
library(nloptr)
library(qs)
library(rlang)
library(readxl)
library(cowplot)

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Documents/ncov_age/covidm/";
cm_force_rebuild = F;
cm_verbose = F;
cm_version = 1;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# Select fIa set to be plotted
fIa = 0.5

# load beijing data
cases_B = fread(path("../Beijing_Onset_20200229.csv"));
cases_B = cases_B[`Date No.` != "Total", .(date = dmy(Date), incidence = as.numeric(Adjusted))];
cases_B = cases_B[date <= ymd("2020-02-22")]

age_dist_B = data.table(read_excel(path("../Beijing_Confirmed_AgeDist.xlsx"), col_names = c("age", "a", "b", "c", "d", "e", "f", "feb22"), skip = 1))[, .(date = "2020-01-01 - 2020-02-22", age, n = feb22)];
age_dist_B[, age := factor(age, levels = unique(age))];
cib = cbind(age_dist_B[, rdirichlet(1000, n + 1) * sum(n), by = date][, .(date, x = V1)], age_dist_B[, rep(age, each = 1000), by = date][, .(age = V1)]);
cib = cib[, cm_mean_hdi(x), keyby = .(date, age)]
age_dist_B = merge(age_dist_B, cib, by = c("date", "age"))
age_dist_B = age_dist_B[, .(age = age, mean = n/sum(n), lower = lower/sum(n), upper = upper / sum(n)), by = date]

# load shanghai data
cases_S = fread(path("../Shanghai_Onset.csv"));
cases_S = cases_S[`Date No.` != "Total", .(date = dmy(Date), incidence = as.numeric(Adjusted))];

age_dist_S = data.table(read_excel(path("../Shanghai_Confirmed_AgeDist.xlsx"), col_names = c("age", "total"), skip = 1, n_max = 2))[, .(date = "2020-01-01 - 2020-02-25", age, n = total)];
age_dist_S[, age := factor(age, levels = unique(age))];
cis = cbind(age_dist_S[, rdirichlet(1000, n + 1) * sum(n), by = date][, .(date, x = V1)], age_dist_S[, rep(age, each = 1000), by = date][, .(age = V1)]);
cis = cis[, cm_mean_hdi(x), keyby = .(date, age)]
age_dist_S = merge(age_dist_S, cis, by = c("date", "age"))
age_dist_S = age_dist_S[, .(age = age, mean = n/sum(n), lower = lower/sum(n), upper = upper / sum(n)), by = date]

# load SK data
patient = fread(path("../SKdataset/patient.csv"))
sk_confirmations = patient[, .(n = .N), keyby = .(date = ymd(confirmed_date))];
sk_confirmations[, population := factor(0)]

sk_age_dist = data.table(age_lower = seq(0, 70, by = 10), age_upper = c(seq(10, 70, by = 10), 100))
sk_age_dist = merge(sk_age_dist, patient[!is.na(birth_year), .N, keyby = .(age_lower = pmin(70, ((2020 - birth_year) %/% 10) * 10))], by = "age_lower");
sk_age_dist[, age := paste(age_lower, "-", age_upper)]
sk_age_dist[, age := factor(age, levels = unique(age))];
sk_age_dist[, date := "all"]
ci = cbind(sk_age_dist[, rdirichlet(1000, N + 1) * sum(N), by = date][, .(date, x = V1)], sk_age_dist[, rep(age, each = 1000), by = date][, .(age = V1)]);
ci = ci[, cm_mean_hdi(x), keyby = .(date, age)]
sk_age_dist = merge(sk_age_dist, ci, by = c("date", "age"))
sk_age_dist = sk_age_dist[, .(age = age, mean = N/sum(N), lower = lower/sum(N), upper = upper / sum(N)), by = date]

# load IT data
matrices = readRDS(path("../all_matrices.rds"));
it_onsets = fread(path("../lombardy-onsets.csv"));
it_confirmations = fread(path("../lombardy-confirmations.csv"));
it_onsets[, date := days_since_feb10 + ymd("2020-02-10")];
it_confirmations[, date := days_since_feb10 + ymd("2020-02-10")];

it_age_dist = data.table(age_lower = seq(0, 70, by = 10), age_upper = c(seq(10, 70, by = 10), 100), N = c(40, 64, 305, 582, 1057, 1760, 1800, 2260 + 1664 + 265));
it_age_dist[, frac := N * 1/sum(N)]
it_age_dist[, age := paste(age_lower, "-", age_upper)]
it_age_dist[, age := factor(age, levels = unique(age))];
it_age_dist[, date := "all"]
ci = cbind(it_age_dist[, rdirichlet(1000, N + 1) * sum(N), by = date][, .(date, x = V1)], it_age_dist[, rep(age, each = 1000), by = date][, .(age = V1)]);
ci = ci[, cm_mean_hdi(x), keyby = .(date, age)]
it_age_dist = merge(it_age_dist, ci, by = c("date", "age"))
it_age_dist = it_age_dist[, .(age = age, mean = N/sum(N), lower = lower/sum(N), upper = upper / sum(N)), by = date]

# load validation fits
bs = cm_load(path(paste0("2-linelist-validation-both-bs-fIa-", fIa, "-con.qs")));
sk = cm_load(path(paste0("2-linelist-validation-both-sk-fIa-", fIa, "-con.qs")));
it = cm_load(path(paste0("2-linelist-validation-both-lb-fIa-", fIa, "-con.qs")));

summ_fit = function(f)
{
    f = melt(f$posterior[, .SD, .SDcols = 5:length(f$posterior)]);
    f = f[, cm_mean_hdi(value), by = variable]
    f = f[, paste0(signif(mean, 2), " (", signif(lower, 2), "-", signif(upper, 2), ")\n")]
    for (i in 1:length(f)) {
        cat(f[i])
    }
}
summ_fit(bs)
summ_fit(sk)
summ_fit(it)

# load posteriors from fitting
fitted = qread(path(paste0("2-linelist_both_fit_fIa", fIa, "-rbzvih.qs")));
fitted_symp_mean = unname(unlist(fitted[, lapply(.SD, mean), .SDcols = y_00:y_70]))
fitted_symp_mean = rep(fitted_symp_mean, each = 2)
fitted_susc_mean = unname(unlist(fitted[, lapply(.SD, mean), .SDcols = u_00:u_70]))
fitted_susc_mean = rep(fitted_susc_mean, each = 2)

# create table
fitted[, .(mean = lapply(.SD, function(x) round(mean(x), 2))), .SDcols = u_00:y_70]
fitted[, .(lapply(.SD, function(x) round(quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)), 2))), .SDcols = u_00:y_70]

bs_d = cm_sample_fit(bs, 100);
sk_d = cm_sample_fit(sk, 100);
it_d = cm_sample_fit(it, 100);

#
bsi = cm_load(path(paste0("2-linelist-validation-both-bs-fIa-", fIa, "-ind.qs")));
ski = cm_load(path(paste0("2-linelist-validation-both-sk-fIa-", fIa, "-ind.qs")));
iti = cm_load(path(paste0("2-linelist-validation-both-lb-fIa-", fIa, "-ind.qs")));

fitted_ind = qread(path(paste0("2-linelist_mn_both_fit_ind_fIa", fIa, "-r.qs")));
fitted_symp_mean = unname(unlist(fitted_ind[location %like% "CCDC", lapply(.SD, mean), .SDcols = y_00:y_70]))
fitted_symp_mean = rep(fitted_symp_mean, each = 2)
fitted_susc_mean = unname(unlist(fitted_ind[location %like% "CCDC", lapply(.SD, mean), .SDcols = u_00:u_70]))
fitted_susc_mean = rep(fitted_susc_mean, each = 2)
bs_di = cm_sample_fit(bsi, 100);

fitted_symp_mean = unname(unlist(fitted_ind[location == "South Korea", lapply(.SD, mean), .SDcols = y_00:y_70]))
fitted_symp_mean = rep(fitted_symp_mean, each = 2)
fitted_susc_mean = unname(unlist(fitted_ind[location == "South Korea", lapply(.SD, mean), .SDcols = u_00:u_70]))
fitted_susc_mean = rep(fitted_susc_mean, each = 2)
sk_di = cm_sample_fit(ski, 100);

fitted_symp_mean = unname(unlist(fitted_ind[location == "Lombardia", lapply(.SD, mean), .SDcols = y_00:y_70]))
fitted_symp_mean = rep(fitted_symp_mean, each = 2)
fitted_susc_mean = unname(unlist(fitted_ind[location == "Lombardia", lapply(.SD, mean), .SDcols = u_00:u_70]))
fitted_susc_mean = rep(fitted_susc_mean, each = 2)
it_di = cm_sample_fit(iti, 100);

# plotting
ll_regions = c("Anhui", "Guangdong", "Guangxi", "Hubei", "Hunan", 
    "Jiangsu", "Jiangxi", "Jilin", "Shaanxi", "Shandong", "Sichuan", 
    "Tianjin", "Zhejiang", "Wuhan_CCDC", "Hubei_CCDC", "ChinaNonHubei_CCDC",
    "Lombardia", "Piemonte", "Trento", "Veneto", 
    "Friuli Venezia Giulia", "Liguria", "Emilia-Romagna", "Toscana", "Marche", 
    "Lazio", "Campania", "Puglia", "Japan", "Singapore", "South Korea", "Ontario")

chin13 = c("#dbdaeb","#c9cfe5","#b4c4df","#9bb9d9","#81aed2","#63a2cb","#4394c3","#2785bb","#1175b0","#0667a1","#045a8d","#034a74","#023858")
ccdc3 = c("#bdd7e7","#6baed6","#2171b5")
ital12 = c("#fdd09c","#fdc28c","#fdb07a","#fb9a66","#f88356","#f16c49","#e65339","#d93826","#c81d13","#b30a06","#9a0101","#7f0000")
singles = c("#7570b3","#e7298a","#66a61e","#e6ab02")

ll_colours = c(chin13, ccdc3, ital12, singles);
names(ll_colours) = ll_regions;

ll_sets = c("Overall", "China (SHO)", "China (CCDC)", "Italy", "Lombardia", "Japan", "Singapore", "South Korea", "Canada")
ll_set_colours = c("#666666", "#4394c3", "#6baed6", "#e65339", "#fdd09c", singles)
names(ll_set_colours) = ll_sets;



modelled = function(dynamics, comp, start_date)
{
    d = dynamics[compartment == comp, .(cases = sum(value)), by = .(t, run, population, scenario)];
    d = d[, cm_mean_hdi(cases), by = .(t, population = as.factor(population), scenario)]
    d[, date := ymd(start_date) + t]
    
    d[, scenario := factor(scenario, ll_sets)]

    ggplot() +    
        geom_ribbon(data = d, aes(date, ymin = lower, ymax = upper, fill = scenario), alpha = 0.25) +
        geom_line(data = d, aes(date, y = mean, colour = scenario), size = 0.25) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = ll_set_colours, aesthetics = c("fill", "colour")) +
        theme(legend.position = "none")
}

add_data = function(data, column, clr)
{
    geom_point(data = data, aes_string(x = "date", y = column), size = 0.2)
}

redo_age = function(x, col = "age")
{
    ages = as.numeric(str_split(x[[col]], " |-|\\+", simplify = T)[,1]);
    ages_lower = unique(ages);
    ages_upper = c(paste0("-", tail(ages_lower, -1) - 1), "+")
    age_ids = match(ages, ages_lower);
    x[, (col) := paste0(..ages_lower[..age_ids], ..ages_upper[..age_ids])];
    x[, (col) := factor(get(col), levels = unique(get(col)))];
    return (x)
}

case_dist = function(dynamics, data, sim_start_date, measurement_end_date, age_bounds, compart = "cases_reported")
{
    res = NULL;
    for (s in dynamics[, unique(scenario)])
    {
        for (i in dynamics[, unique(run)])
        {
            cat(".");
            res = rbind(res,
                cbind(cm_case_distribution(dynamics[scenario == s & run == i],
                    sim_start_date, sim_start_date, measurement_end_date, age_bounds, compart), scenario = s, run = i));
        }
    }
    cat("\n");

    data = rbind(
        cbind(res[, cm_mean_hdi(fcases), by = .(age, scenario)]),
        cbind(data[, .(age, mean, lower, upper)], scenario = "observed")
    );
    
    data[, scenario := factor(scenario, c("observed", ll_sets))];
    data = redo_age(data);
    
    ggplot(data) +
        geom_col(data = data[scenario != "observed"], 
            aes(x = age, y = 100 * upper, fill = scenario, group = scenario), colour = NA, alpha = 0.5,
            position = position_dodge(width = 0.8), width = 0.7) +
        geom_col(data = data[scenario != "observed"], 
            aes(x = age, y = 100 * lower, fill = scenario, group = scenario), colour = NA, 
            position = position_dodge(width = 0.8), width = 0.7) +
        geom_point(data = data[scenario != "observed"], 
            aes(x = age, y = 100 * mean, colour = scenario, group = scenario), 
            size = 0.25, position = position_dodge(width = 0.8), shape = 20) +
        geom_col(data = data[scenario == "observed"], 
            aes(x = age, y = 100 * mean, group = scenario), fill = NA, colour = "black", 
            size = 0.125, position = position_dodge(width = 0.9)) +
        geom_errorbar(data = data[scenario == "observed"], 
            aes(x = age, ymin = 100 * lower, ymax = 100 * upper, group = scenario), 
            size = 0.125, width = 0.2, position = position_dodge(width = 0.9), colour = "black") +
        scale_fill_manual(values = ll_set_colours, aesthetics = c("fill", "colour")) +
        labs(x = "Age", y = "Clinical cases (%)") +
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plotting setup
txt_theme = theme(plot.title = element_text(face = "plain", size = 7, hjust = 0.5))
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# Stitch together validation fits
bs_all = rbind(
    cbind(bs_d, scenario = "Overall"),
    cbind(bs_di, scenario = "China (CCDC)")
)

it_all = rbind(
    cbind(it_d, scenario = "Overall"),
    cbind(it_di, scenario = "Lombardia")
)

sk_all = rbind(
    cbind(sk_d, scenario = "Overall"),
    cbind(sk_di, scenario = "South Korea")
)

epi1 = modelled(bs_all[population == 1], "cases_reported", bs$base_parameters$date0) +
    add_data(cases_B, "incidence") + labs(x = "Date", y = "Incident cases", title = "Beijing") + txt_theme

epi2 = modelled(bs_all[population == 2], "cases_reported", bs$base_parameters$date0) +
    add_data(cases_S, "incidence") + labs(x = "Date", y = NULL, title = "Shanghai") + txt_theme

epi3 = modelled(sk_all[t > 15], "cases_reported", sk$base_parameters$date0) +
    add_data(sk_confirmations, "n") + labs(x = "Date", y = NULL, title = "South Korea") + txt_theme

epi4 = modelled(it_all[t > 30], "cases_reported", it$base_parameters$date0) +
    add_data(it_confirmations, "confirmations") + labs(x = "Date", y = NULL, title = "Lombardy") + txt_theme

dist1 = case_dist(bs_all[population == 1], age_dist_B, bs$base_parameters$date0, "2020-03-01", c(0, 6, 18, 60, 100))
dist2 = case_dist(bs_all[population == 2], age_dist_S, bs$base_parameters$date0, "2020-03-01", c(0, 18, 100)) + labs(y = NULL)
dist3 = case_dist(sk_all, sk_age_dist, sk$base_parameters$date0, "2020-03-15", c(seq(0, 70, by = 10), 100)) + labs(y = NULL)
dist4 = case_dist(it_all, it_age_dist, it$base_parameters$date0, "2020-03-15", c(seq(0, 70, by = 10), 100)) + labs(y = NULL)

f2c = plot_grid(epi1, epi2, epi3, epi4,
          dist1, dist2, dist3, dist4,
          nrow = 2, ncol = 4, rel_heights = c(1, 1), rel_widths = c(9, 8, 8, 8))
f2c
ggsave(path(paste0("../Submission/Fig2c_fIa", fIa, ".pdf")), f2c, width = 10, height = 8, units = "cm", useDingbats = F)


# Show constituent age distributions
get_ind = function(fIa)
{
  age_dist = fread(path("age.csv"));
  ind = qread(path(paste0("2-linelist_mn_both_fit_ind_fIa", fIa, "-r.qs")));
  ind_expected = qread(path(paste0("2-linelist_both_expected_fIa", fIa, "-r.qs")));
  
  ind_expected = melt(ind_expected, id.vars = "location")
  ind_expected[, age_lower := as.numeric(str_sub(variable, 3))]
  ind_expected[, age_upper := age_lower + 9]
  ind_expected[age_upper == 79, age_upper := 119]
  ind_expected[, age := paste(age_lower, "-", age_upper)]
  ind_expected = ind_expected[, cm_mean_hdi(value), by = .(location, age)]
  
  # amalgamate over 70s
  age_dist = rbind(age_dist[age_lower < 70], age_dist[age_lower >= 70, .(age_lower = 70, age_upper = 120, n = sum(n)), by = location])
  
  age_dist[, age := paste(age_lower, "-", age_upper - 1)]
  age_dist[, age := factor(age, levels = unique(age))];
  ci = cbind(age_dist[, rdirichlet(1000, n + 1) * sum(n), by = location][, .(location, x = V1)], age_dist[, rep(age, each = 1000), by = location][, .(age = V1)]);
  ci = ci[, cm_mean_hdi(x), keyby = .(location, age)]
  age_dist = merge(age_dist, ci, by = c("location", "age"))
  age_dist = age_dist[, .(age = age, mean = n/sum(n), lower = lower/sum(n), upper = upper / sum(n)), by = location]
  
  ind = merge(age_dist[, .(location, age, cmean = mean, clower = lower, cupper = upper)],
        ind_expected[, .(location, age, emean = mean, elower = lower, eupper = upper)],
        by = c("location", "age"))
  
  ind[, age := factor(age, levels = unique(age))]
  ind[, location := factor(location, levels = ll_regions)]
  ind[, fIa := fIa];
  
  return (ind)
}

theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

ind = get_ind(fIa)
ind = redo_age(ind)
f2a = ggplot(ind) +
    geom_col(aes(x = age, y = cmean), fill = "#ffffff", colour = "#000000", size = 0.125) +
    geom_errorbar(aes(x = age, ymin = clower, ymax = cupper), size = 0.125, width = 0.2) +
    geom_ribbon(aes(as.numeric(age), ymin = elower, ymax = eupper, fill = location), alpha = 0.5) +
    geom_line(aes(as.numeric(age), y = emean, colour = location)) +
    facet_wrap(~location, ncol = 4) + 
    labs(x = "Age", y = "Reported cases") +
    scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.5)) +
    scale_fill_manual(values = ll_colours, aesthetics = c("fill", "colour")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), legend.position = "none")

ggsave(paste0("~/Dropbox/nCoV/Submission/Fig2a_fIa", fIa, ".pdf"), f2a, width = 10, height = 14, unit = "cm", useDingbats = F)


# symp rate
ind = qread(path(paste0("2-linelist_mn_both_fit_ind_fIa", fIa, "-r.qs")));
ind[, location := factor(location, levels = ll_regions)]

ind = melt(ind[, c(5:20, 23)], id.vars = "location")
ind[, var := str_sub(variable, 1, 1)];
ind[, age := as.numeric(str_sub(variable, 3))]
ind[, variable := NULL]
ind2 = ind[, cm_mean_hdi(value), by = .(location, var, age)]
fS = ggplot(ind2) + 
    geom_ribbon(aes(x = age + 5, ymin = lower, ymax = upper, fill = location, group = var), alpha = 0.25) +
    geom_line(aes(x = age + 5, y = mean, colour = location, linetype = var), size = 0.25) +
    labs(x = "Age", y = "Clinical fraction") +
    facet_wrap(~location, ncol = 4) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_linetype_manual(values = c("u" = "dashed", "y" = "solid")) +
    scale_fill_manual(values = ll_colours, aesthetics = c("fill", "colour")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), legend.position = "none")

ggsave(paste0("~/Dropbox/nCoV/Submission/FigS_symp_fIa", fIa, ".pdf"), fS, width = 10, height = 14, unit = "cm", useDingbats = F)


# overall - run symp rate first
isets = data.table(location = ll_regions, set = c(rep("China (SHO)", 13), rep("China (CCDC)", 3), rep("Italy", 12), "Japan", "Singapore", "South Korea", "Canada"))
ind3 = merge(ind, isets, by = "location");
ind3 = ind3[, cm_mean_hdi(value), by = .(set, var, age)]

fitted_symp = qread(path(paste0("2-linelist_both_fit_fIa", fIa, "-rbzvih.qs")));
fs = melt(fitted_symp[, 5:20], id.vars = NULL)
fs[, var := str_sub(variable, 1, 1)];
fs[, age := as.numeric(str_sub(variable, 3))]
fs2 = fs[, cm_mean_hdi(value), by = .(var, age)]

s_all = rbind(ind3, cbind(set = "Overall", fs2))
s_all[, set := factor(set, levels = ll_sets)]
s_all[var == "u", var := "Susceptibility"]
s_all[var == "y", var := "Clinical fraction"]

f2b = ggplot(s_all) + 
    geom_ribbon(aes(x = age + 5, ymin = lower, ymax = upper, fill = set, group = var), alpha = 0.25) +
    geom_line(aes(x = age + 5, y = mean, colour = set, linetype = var), size = 0.25) +
    labs(x = "Age", y = "Clinical fraction + susceptibility", linetype = NULL) +
    facet_wrap(~set, nrow = 2, ncol = 4) +
    scale_linetype_manual(values = c("Susceptibility" = "dotted", "Clinical fraction" = "solid")) +
    scale_fill_manual(values = ll_set_colours, aesthetics = c("fill", "colour")) +
    theme(strip.background = element_blank(), legend.position = c(0.02, 0.6)) +
    ylim(0, NA) + xlim(0, 80) +
    guides(fill = "none", colour = "none", linetype = "legend")
ggsave(paste0("~/Dropbox/nCoV/Submission/Fig2b_fIa", fIa, ".pdf"), f2b, width = 12, height = 6, unit = "cm", useDingbats = F)

# overall figure 2
f2 = plot_grid(f2a, plot_grid(f2b, f2c, ncol = 1, rel_heights = c(6, 8), labels = c("b", "c"), label_size = 10), nrow = 1, rel_widths = c(10, 12, 8), labels = c("a", "", ""), label_size = 10)
ggsave(paste0("~/Dropbox/nCoV/Submission/Fig2_fIa", fIa, ".pdf"), f2, width = 22, height = 14, unit = "cm", useDingbats = F)




# POSTERIORS
p1 = cm_plot_posterior(cm_load(path("2-linelist-validation-both-bs-fIa-0.5-ind.qs")))
p2 = cm_plot_posterior(cm_load(path("2-linelist-validation-both-sk-fIa-0.5-ind.qs")))
p3 = cm_plot_posterior(cm_load(path("2-linelist-validation-both-lb-fIa-0.5-ind.qs")))
p = plot_grid(p1, p2, p3, ncol = 1, labels = c("a", "b", "c"), label_size = 8, rel_heights = c(3, 2.2, 3))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/S-2-validation-posteriors.pdf", p, width = 12, height = 15, units = "cm", useDingbats = F)

q = cm_load(path("2-linelist_both_fit_fIa0.5-rbzvih.qs"))
q = melt(q, measure.vars = 5:23, id.vars = "ll")
p = ggplot(q) + geom_histogram(aes(x = value, fill = variable), bins = 30) + facet_wrap(~variable, scales = "free", ncol = 6) + 
    theme(legend.position = "none", strip.background = element_blank())
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/S-posterior-both-fit-con.pdf", p, width = 20, height = 8, units = "cm", useDingbats = F)




# how clinical fraction changes with fIa plot

data = NULL
for (fIa in c(0, 0.25, 0.5, 0.75, 1.0)) {
    fitted_symp = qread(path(paste0("2-linelist_both_fit_fIa", fIa, "-rbzvih.qs")));
    fs = melt(fitted_symp[, 5:20], id.vars = NULL)
    fs[, age := as.numeric(str_sub(variable, 3))]
    fs[, var := str_sub(variable, 1, 1)]
    fs[, fIa := fIa]
    fs2 = fs[, cm_mean_hdi(value), by = .(age, var, fIa)]

    data = rbind(data, cbind(set = "Overall", fs2))
}

data[, set := factor(set, levels = ll_sets)]
data[var == "u", var := "Susceptibility"]
data[var == "y", var := "Clinical fraction"]

big_fig = ggplot(data) + 
    geom_line(aes(x = age + 5, y = mean, linetype = as.factor(fIa), group = fIa), size = 0.3) +
    geom_ribbon(aes(x = age + 5, ymin = lower, ymax = upper, group = fIa), alpha = 0.125) +
    labs(x = "Age", y = "Clinical fraction", linetype = "Subclinical infectiousness") +
    facet_wrap(~var, nrow = 2, ncol = 4) +
    scale_linetype_manual(values = c("solid", "82", "4212", "22", "11")) +
    scale_fill_manual(values = ll_set_colours, aesthetics = c("fill", "colour")) +
    theme(strip.background = element_blank(), legend.key.width = unit(1, "cm"), legend.position = "bottom") +
    ylim(0, 1)

ggsave(paste0("~/Dropbox/nCoV/Submission/Supp Figs/FigS-fIa_symp.pdf"), big_fig, width = 12, height = 5, unit = "cm", useDingbats = F)

