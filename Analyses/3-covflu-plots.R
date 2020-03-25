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
library(spatstat)

txt_theme = theme(plot.title = element_text(face = "plain", size = 7, hjust = 0.5))
path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

ccols = c("Milan" = "#623b84", "Birmingham" = "#518c9c", "Bulawayo" = "#a3d55a")
c_colours = function(a = c("colour", "fill")) scale_colour_manual(values = ccols, aesthetics = a)


# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
cm_verbose = F
if (Sys.info()["nodename"] %like% "lshtm") {
   cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# Return the mean and CI from a set of samples x.
mean_ci = function(x, ci = 0.95)
{
    h = quantile(x, c(0.5 - ci/2, 0.5 + ci/2));
    m = mean(x);
    return (list(mean = m, lower = h[[1]], upper = h[[2]]))
}

# Create flu scenario plot
flu_scenario = data.table(read_excel(path("../flu_scenario_v3.xlsx")));
flu_scenario[, Susceptibility := OldSusc]
flu_scenario[, Severity := scaled_newSev]
flu_scenario[, Age_band := factor(Age_band, levels = unique(Age_band))]

# Load and process data
argv = commandArgs();
argc = length(argv);
a_which = argv[argc];

covflu = cm_load(path(paste0("3-covflu-fix-", a_which, ".qs")))
fwrite(covflu$R0[, mean(R0), keyby = .(location, virus, fIa)], paste0("~/Dropbox/nCoV/Analyses/3-covflu-", a_which, "-R0.csv"));

#epi = results[, cm_mean_hdi(cases), by = .(t, virus, location, fIa, schools)]
epi = covflu$epi[, cm_mean_hdi(cases), by = .(t, virus, location, fIa, schools, popK)]
epi[, schools := factor(schools, levels = c("open", "closed"))]

dist = covflu$dist
#Plotting
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# Plot a population pyramid
popdist = function(populations, countries)
{
    pop = populations[name %in% countries];
    pop2 = NULL
    
    for (co in countries)
    {
        p = pop[name == co];
        f85p = p[name == co & age %like% "85|90|95|100", sum(f)]
        m85p = p[name == co & age %like% "85|90|95|100", sum(m)]
        p = p[! (age %like% "85|90|95|100")];
        p = rbind(p, p[.N]);
        p[.N, age := "85+"]
        p[.N, f := f85p]
        p[.N, m := m85p]
        pop2 = rbind(pop2, p);
    }
    
    pop2[, name := str_replace(name, "[A-Za-z]+ \\| ", "")]
    
    yl = "Population (thousands)";

    median_age = pop2[, weighted.median(1:.N, f+m) * 5 - 2.5, by = name]

    median_age[, name := factor(name, levels = unique(name))]
    pop2[, name := factor(name, levels = unique(name))]
    
    ages = unique(pop2$age);

    ggplot(pop2) +
        geom_col(aes(x = age, y = f + m, fill = name), alpha = 0.5) +
        geom_text(data = median_age, aes(label = paste("Median age:", round(V1))), x = 11, y = 10, hjust = 0, size = 2) +
        coord_flip() +
        scale_x_discrete(limits = ages, breaks = ages[seq(1, length(ages), by = 2)]) +
        facet_wrap(~name, nrow = 1) +
        labs(x = NULL, y = yl) +
        theme(strip.background = element_blank(), legend.position = "none") +
        c_colours() +
        txt_theme
}

# version that shows only one fIa
epicol = function(covflu, epi, vir, Ia, titl) {
    epi2 = covflu$epi[, .(total = sum(cases)), by = .(location, virus, schools, fIa, run, popK)]
    epi2 = epi2[, cm_mean_hdi(total), by = .(location, virus, schools, fIa, popK)]
    epi2[, text := paste0(schools, ": ", round(mean/1000), "K (", round(lower/1000), "-", round(upper/1000), ")")]
    th = epi[fIa == Ia & virus == vir, max(upper/popK)]
    epi2[, textheight := th]
    epi2[schools == "closed", textheight := textheight * 0.88];
    
    ggplot(epi[fIa == Ia & virus == vir]) + 
        geom_ribbon(aes(x = t, ymin = lower/popK, ymax = upper/popK, group = schools), alpha = 0.25) +
        geom_line(aes(x = t, y = mean/popK, group = schools, linetype = schools)) +
        geom_text(data = epi2[fIa == Ia & virus == vir], aes(x = 225, y = textheight, label = text), size = 1.5) +
        facet_wrap(~location, ncol = 1) +
        c_colours() +
        labs(x = "Time (days)", y = "Clinical cases per 1000 population", title = titl) +
        theme(legend.position = c(0.5, 0.82))
}

# a DIAGRAM
covid_scenario = qread(path(paste0("2-linelist_symp_fit_fIa0.5.qs")));
covid_scenario = colMeans(covid_scenario[, 5:12])
covid_scenario = as.data.table(as.list(covid_scenario))
covid_scenario = melt(covid_scenario, variable.name = "age", value.name = "Severity")
covid_scenario[, age := as.numeric(str_sub(age, 3)) + 5]
covid_scenario[, Susceptibility := 0.5]
covid_scenario = rbind(covid_scenario, covid_scenario[8, .(age = 95, Severity, Susceptibility)])
flu_scenario[, age := Age_gp * 5 - 2.5]

scen = rbind(
    covid_scenario[, .(age, Severity, Susceptibility, virus = "COVID")],
    flu_scenario[, .(age, Severity, Susceptibility, virus = "Influenza-like")]
)

f3a = ggplot(scen) +
    geom_line(aes(x = age, y = Severity, group = 1, colour = "Severity")) +
    geom_line(aes(x = age, y = Susceptibility, group = 2, colour = "Susceptibility")) +
    facet_wrap(~virus, nrow = 1) +
    labs(x = "Age", y = NULL, colour = NULL) + ylim(0,1) + xlim(0, 100) +
    theme(legend.position = c(0.025, 0.88), strip.background = element_blank()) + txt_theme

ggsave(path(paste0("../Submission/Fig3a-", a_which, ".pdf")), f3a, width = 7, height = 3, units = "cm", useDingbats = F)


# b POP PYRAMIDS
f3b = popdist(cm_populations, c("Italy | Milan", "UK | Birmingham", "Zimbabwe | Bulawayo"))
ggsave(path(paste0("../Submission/Fig3b-", a_which, ".pdf")), f3b, width = 7, height = 3.5, units = "cm", useDingbats = F)

# c ATTACK RATE BY AGE
popgs = NULL;
for (nm in c("Italy | Milan", "UK | Birmingham", "Zimbabwe | Bulawayo")) {
    popcat = cm_populations[name == nm]
    popcat[, age_lower := as.numeric(str_replace_all(age, "(\\+|-[0-9]+)", ""))]
    popcat[, age_group := age_lower %/% 10]
    popcat[age_group >= 7, age_group := 7]
    popgs = rbind(popgs, popcat[, .(location = nm, popKg = sum(f+m)), by = age_group]);
}
popgs[, location := str_replace_all(location, "[^|]*\\| ", "")]
popgs[, age := paste(age_group * 10, "-", (age_group + 1) * 10)]
popgs[, age_group := NULL]

dist[, schools := factor(schools, levels = c("open", "closed"))]
dist = merge(dist, unique(epi[, .(location, popK)]), by = c("location"))
dist = merge(dist, popgs, by = c("location", "age"))
dist[, location := factor(location, levels = c("Milan", "Birmingham", "Bulawayo"))]

dist2 = dist[, cm_mean_hdi(cases/popK), by = .(age, virus, location, fIa, schools)];

f3c1 = ggplot(dist2[fIa == 0.5 & schools == "open"]) + 
    geom_col(aes(age, mean, fill = location)) +
    geom_errorbar(aes(age, ymin = lower, ymax = upper), width = 0.2, size = 0.125) +
    facet_grid(virus~location, scales = "free_y") +
    labs(x = "Age", y = "Attack rate\nper 1000 population") +
    c_colours() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), strip.text = element_blank())

f3c = ggdraw(f3c1) + 
    draw_text("COVID", 0.18, 0.95, size = 7, hjust = 0) +
    draw_text("Influenza-like", 0.18, 0.57, size = 7, hjust = 0)

ggsave(path(paste0("../Submission/Fig3c-", a_which, ".pdf")), f3c, width = 7, height = 5, units = "cm", useDingbats = F)


# d EPIDEMICS
# plot_grid(epicol(covflu, epi, "cov", 0.5, "COVID"), 
#           epicol(covflu, epi, "flu", 0.5, "Influenza-like"), nrow = 1)
# 

# version that splits by fIa
epicol2 = function(epi, vir, Ia, titl, yl, lp = "none") {
    ggplot(epi[fIa %in% Ia & virus == vir & t < 250]) + 
        geom_ribbon(aes(x = t, ymin = lower/popK, ymax = upper/popK, group = paste(schools,fIa), fill = location), alpha = 0.25) +
        geom_line(aes(x = t, y = mean/popK, group = paste(schools,fIa), linetype = schools, colour = location), size = 0.25) +
        facet_grid(fIa ~ location) +
        labs(x = "Time (days)", y = yl, title = titl) +
        guides(colour = F, fill = F) +
        c_colours() +
        theme(legend.position = lp, strip.background = element_blank(), 
            legend.text = element_text(size = 6), legend.title = element_text(size = 6), 
            legend.spacing.y = unit(0.05, "cm"), legend.key.height = unit(0.2, "cm")) + txt_theme +
        scale_x_continuous(breaks = c(0, 100, 200))
}

f3d = plot_grid(epicol2(epi, "cov", c(0.0, 0.25, 0.5, 0.75), "COVID", "Clinical cases\nper 1000 population", c(0.8, 0.92)),
    epicol2(epi, "flu", c(0.0, 0.25, 0.5, 0.75), "Influenza-like", NULL), nrow = 1)

ggsave(path(paste0("../Submission/Fig3d-", a_which, ".pdf")), f3d, width = 10, height = 5, units = "cm", useDingbats = F)


# e IMPACT OF SCHOOL CLOSURES
stats = covflu$epi[, .(peak_time = t[which.max(cases)], peak_cases = max(cases), total_cases = sum(cases)), 
                   by = .(run, virus, location, fIa, schools, popK)]

stats2 = stats[, .(peak_time = median(peak_time), peak_cases = median(peak_cases), total_cases = median(total_cases)), 
    by = .(virus, location, fIa, schools, popK)]
stats2 = cbind(stats2[schools == "open"], stats2[schools == "closed"])
names(stats2)[9:16] = paste0("cl", names(stats2)[9:16])

stats[, peak_time := peak_time + rnorm(.N, 0, 0.25)]
stats[virus == "flu", virus := "Influenza-like"]
stats[virus == "cov", virus := "COVID"]
stats2[virus == "flu", virus := "Influenza-like"]
stats2[virus == "cov", virus := "COVID"]

f3e = ggplot(stats[fIa == 0.5], aes(peak_time, peak_cases / popK)) +
    geom_segment(data = stats2[fIa == 0.5],
             aes(x = peak_time, y = peak_cases / popK, xend = clpeak_time, yend = clpeak_cases / clpopK), colour = "black", size = 0.125) +
    geom_point(aes(colour = location), alpha = 0.5, size = 0.5) +
    stat_ellipse(aes(linetype = schools, group = paste(location, schools, virus)), colour = "black", size = 0.25, type = "norm", alpha = 0.75, level = 0.97) +
    labs(x = "Time to peak (days)", y = "Peak incidence\nper 1000 population") +
    ylim(0, NA) +
    c_colours() + guides(colour = FALSE) +
    theme(legend.position = c(0.68, 0.7)) +
    facet_wrap(~virus, scales = "free_y", nrow = 2) +
    theme(strip.background = element_blank(), legend.text = element_text(size = 6), legend.title = element_text(size = 6), 
            legend.spacing.y = unit(0.05, "cm"), legend.key.height = unit(0.2, "cm")) + txt_theme

ggsave(path(paste0("../Submission/Fig3e-", a_which, ".pdf")), f3e, width = 10, height = 6, units = "cm", useDingbats = F)

# # total cases version
# f3d = ggplot(stats[fIa == 0.5], aes(peak_time, total_cases / popK)) +
#     geom_point(aes(colour = location), alpha = 0.1, size = 0.5) +
#     stat_ellipse(aes(colour = location, linetype = schools, group = paste(location, schools, virus)), type = "norm", alpha = 0.75) +
#     labs(x = "Time to peak", y = "Total cases (per 1000 people)") +
#     theme(legend.position = c(0.8, 0.6)) +
#     geom_segment(data = stats2[fIa == 0.5],
#              aes(x = peak_time, y = total_cases / popK, xend = clpeak_time, yend = cltotal_cases / clpopK, colour = location),
#                  linetype = "11") + scale_y_log10()
# 



# f IMPACT OF SCHOOL CLOSURES
stats = covflu$epi[, .(peak_time = t[which.max(cases)], peak_cases = max(cases), total_cases = sum(cases)), 
                   by = .(run, virus, location, fIa, schools, popK)]
stats[virus == "flu", virus := "Influenza-like"]
stats[virus == "cov", virus := "COVID"]

stats2 = stats[, .(peak_time = median(peak_time), peak_cases = median(peak_cases), total_cases = median(total_cases)), 
    by = .(virus, location, fIa, schools, popK)]
stats3 = cbind(stats2[schools == "open"], stats2[schools == "closed"])
names(stats3)[9:16] = paste0("cl", names(stats3)[9:16])

dnudge = stats3[virus == "COVID", max(peak_cases / popK)] / 10
if (a_which == "R0") {
    dnudge = dnudge * 0.6
}

f3f = ggplot() +
    geom_segment(data = stats3[virus == "COVID"],
        aes(x = peak_time, y = peak_cases / popK, xend = clpeak_time, yend = clpeak_cases / clpopK), colour = "black", size = 0.125, linetype = "93") +
    geom_text(data = stats3[virus == "COVID" & fIa < 0.6],
        aes(x = (peak_time + clpeak_time)/2, y = (peak_cases / popK + clpeak_cases / clpopK) / 2, label = fIa), size = 1.5, colour = "black", nudge_x = 0.8*dnudge, nudge_y = dnudge/4, hjust = 0) +
    geom_text(data = stats3[virus == "COVID" & fIa > 0.6],
        aes(x = (peak_time + clpeak_time)/2, y = (peak_cases / popK + clpeak_cases / clpopK) / 2, label = fIa), size = 1.5, colour = "black", nudge_x = -0.8*dnudge, nudge_y = -dnudge/4, hjust = 1) +
    geom_point(data = stats2[virus == "COVID"], 
        aes(x = peak_time, y = peak_cases / popK, colour = location, shape = schools), size = 1.5, alpha = 0.75) +
    labs(x = "Time to peak (days)", y = "Peak incidence\nper 1000 population") +
    facet_wrap(~virus, scales = "free_y", nrow = 2) +
    theme(strip.background = element_blank()) + txt_theme + guides(colour = FALSE) +
    c_colours() +
    theme(legend.position = c(0.68, 0.7), strip.background = element_blank(), legend.text = element_text(size = 6), legend.title = element_text(size = 6), 
            legend.spacing.y = unit(0.05, "cm"), legend.key.height = unit(0.2, "cm"))

ggsave(path(paste0("../Submission/Fig3f-", a_which, ".pdf")), f3f, width = 5, height = 6.5, units = "cm", useDingbats = F)


f3 = plot_grid(
    plot_grid(f3a, f3b, f3c, rel_heights = c(3, 3.5, 5), ncol = 1, labels = c("a", "b", "c"), label_size = 9), 
    plot_grid(f3d, 
        plot_grid(f3e, f3f, nrow = 1, labels = c("e", "f"), label_size = 9),
        rel_heights = c(5, 4), ncol = 1, labels = c("d", "", ""), label_size = 9),
    nrow = 1, rel_widths = c(7, 13))
ggsave(path(paste0("../Submission/Fig3-", a_which, ".pdf")), f3, width = 20, height = 10, units = "cm", useDingbats = F)

