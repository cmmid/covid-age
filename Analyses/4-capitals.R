library(ggplot2)
library(data.table)
library(lubridate)
library(stringr)
library(readxl)
library(spatstat)

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_version = 1
cm_force_rebuild = F;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_builddir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

flu_scenario = data.table(read_excel(path("../flu_scenario_v3.xlsx")));
flu_scenario[, Susceptibility := OldSusc]
flu_scenario[, Severity := scaled_newSev]
covid_scenario = qread(path("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
variability = mean(unlist(covid_scenario[, lapply(.SD, function(x) sd(x)/mean(x)), .SDcols = y_00:y_70]));
matrices = readRDS(path("../all_matrices.rds"));

# for setting R0
calc_R0 = function(p, contact, susc_frac = 1) {
    po = p$pop[[1]];
    dIp = sum(po$dIp * seq(0, by = p$time_step, length.out = length(po$dIp)));
    dIs = sum(po$dIs * seq(0, by = p$time_step, length.out = length(po$dIs)));
    dIa = sum(po$dIa * seq(0, by = p$time_step, length.out = length(po$dIa)));
    
    ngm = susc_frac * po$u * t(t(contact) * (
        po$y * (po$fIp * dIp + po$fIs * dIs) + 
        (1 - po$y) * po$fIa * dIa)
    )
    abs(eigen(ngm)$values[1])
}

# parameters
parameters = function(pop_region, mat_region, virus, R0, fIa, report, qS, lower_income)
{
    # Preparing country-specific parameters
    mats = matrices[[mat_region]];
    n_age_groups = nrow(mats$home);
    pop_size = cm_make_population(pop_region, n_age_groups);
    
    # set global parameters
    p = list();
    p$time_step = 0.25;
    p$date0 = "2020-01-01";
    p$time0 = "2020-01-01";
    p$time1 = "2022-01-01";
    p$report_every = 4;
    p$travel = matrix(1.0, nrow = 1, ncol = 1);
    p$fast_multinomial = F;
    p$deterministic = T;
    
    # set population parameters
    pop = list();
    pop$type = "SEI3R";
    
    pop$dE  = cm_delay_gamma(3.0, 4.0, t_max = 60, t_step = 0.25)$p # Derived from He et al (44% infectiousness presymptomatic) https://www.nature.com/articles/s41591-020-0869-5#Sec9
    pop$dIp = cm_delay_gamma(2.1, 4.0, t_max = 60, t_step = 0.25)$p # Derived from Lauer et al (5.1 day incubation period) https://www.ncbi.nlm.nih.gov/pubmed/32150748
    pop$dIs = cm_delay_gamma(2.9, 4.0, t_max = 60, t_step = 0.25)$p # 5 days total: 5.5 day serial interval
    pop$dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p # Assumed same infectious period as clinical cases
    pop$dH  = c(1, 0); # hospitalization ignored
    pop$dC  = c(1, 0); # cases are reported immediately -- because fitting is retrospective to case onset.
    
    pop$size = pop_size;
    pop$matrices = mats;
    pop$contact = c(1, 1, 1, 1);
    pop$contact_mult = c(1, 1, 1, 1);
    pop$contact_lowerto = c(100, 100, 100, 100);
    pop$fIp = rep(1, n_age_groups);
    pop$fIs = rep(1, n_age_groups);
    pop$fIa = rep(fIa, n_age_groups)
    pop$rho = rep(report, n_age_groups);
    pop$tau = rep(1, n_age_groups);
    pop$seed_times = rep(0:27, each = 2);
    pop$dist_seed_ages = cm_age_coefficients(20, 50, 5 * 0:n_age_groups);
    pop$schedule = list()
    
    # TODO with fewer age groups this doesn't quite work nicely
    if (virus == "flu") {
        pop$u = c(flu_scenario$Susceptibility[1:(n_age_groups - 1)], mean(flu_scenario$Susceptibility[n_age_groups:nrow(flu_scenario)]));
        pop$y = c(flu_scenario$Severity[1:(n_age_groups - 1)], mean(flu_scenario$Severity[n_age_groups:nrow(flu_scenario)]));
        #pop$y = pmax(0, pmin(1, pop$y * rnorm(length(pop$y), 1, variability)));
    } else {
        # # MEAN VERSION
        # covy = unname(unlist(covid_scenario[, colMeans(.SD), .SDcols = y_00:y_70]));
        # covu = unname(unlist(covid_scenario[, colMeans(.SD), .SDcols = u_00:u_70]));
        # SAMPLING VERSION
        covy = unname(unlist(covid_scenario[sample.int(nrow(covid_scenario), 1), y_00:y_70]));
        covu = unname(unlist(covid_scenario[sample.int(nrow(covid_scenario), 1), u_00:u_70]));
        
        if (lower_income) {
            covy = covy[c(1,2,4,5,6,7,8,8)];
            covy[1:2] = covy[1:2] + 0.15
        }
        
        covy = unname(rep(covy, each = 2))
        covu = unname(rep(covu, each = 2))
        if (length(covy) < n_age_groups) {
            covy = c(covy, rep(tail(covy, 1), n_age_groups - length(covy)));
            covu = c(covu, rep(tail(covu, 1), n_age_groups - length(covu)));
        } else if (n_age_groups < length(covy)) {
            covy = c(covy[1:(n_age_groups - 1)], mean(covy[n_age_groups:length(covy)]));
            covu = c(covu[1:(n_age_groups - 1)], mean(covu[n_age_groups:length(covu)]));
        }
        pop$y = covy;
        pop$u = covu * 0.1;
    }

    # calculated in 3-covflu.R
    if (fIa == 0) {
        pop$u = 0.1055 * pop$u / mean(pop$u);
    } else if (fIa == 0.5) {
        pop$u = 0.0627 * pop$u / mean(pop$u);
    } else if (fIa == 1.0) {
        pop$u = 0.0429 * pop$u / mean(pop$u);
    } else {
        stop("only fIa in 0, 0.5, 1 supported");
    }
    
    p$pop[[1]] = pop;
    
    return (cm_translate_parameters(p))
}


# Diagram for low-income shift
if (0) {
    yhigh = parameters("Italy", "Italy", "cov", -1, 0.5, 1, 1, F)$pop[[1]]$y;
    ylow = parameters("Italy", "Italy", "cov", -1, 0.5, 1, 1, T)$pop[[1]]$y;

    illustration = data.table(age = rep(seq(5, 80, by = 5), 2), y = c(yhigh, ylow), income = rep(c("high and upper middle", "lower middle and low"), each = 16))
    ggplot(illustration) +
        geom_step(aes(x = age, y = y, colour = income), direction = "vh") +
        labs(x = "Age", y = "Clinical proportion") +
        theme(panel.grid = element_line(colour = "#dddddd")) + ylim(0, 1)
    ggsave("~/Dropbox/nCoV/Submission/Supp Figs/low_income_adjust.pdf", width = 10, height = 3, units = "cm", useDingbats = F)
}

shrink = function(dyn)
{
    setDT(dyn)
    dyn[compartment %in% c("cases", "subclinical"), .(cases = sum(value)), by = .(t, compartment)]
}

cased = function(dyn, run, virus, location, fIa, schools)
{
    setDT(dyn)
    dyn2 = dyn[compartment == "cases", .(total = sum(value)), by = t]
    dyn2[, total := pmax(0, total)] # TODO why does this go negative? deterministic mode may just be using compartment.size which can accumulate rounding errors
    dyn2[, progression := cumsum(total) / sum(total)]
    breaks = findInterval(c(1/3, 2/3), dyn2$progression);
    break1 = ymd("2020-01-01") + breaks[1];
    break2 = ymd("2020-01-01") + breaks[2];

    rbind(
        cbind(cm_case_distribution(dyn, "2020-01-01", "2020-01-01", break1, seq(0, 80, by = 10)),
              run = run, virus = virus, location = location, fIa = fIa, schools = schools, period = "early"),
        cbind(cm_case_distribution(dyn, "2020-01-01", break1, break2, seq(0, 80, by = 10)),
              run = run, virus = virus, location = location, fIa = fIa, schools = schools, period = "middle"),
        cbind(cm_case_distribution(dyn, "2020-01-01", break2, "2021-01-01", seq(0, 80, by = 10)),
              run = run, virus = virus, location = location, fIa = fIa, schools = schools, period = "late")
    )
}



#
# WORLD CAPITALS, DEMOGRAPHICS & MATRICES
#

library(maps)
library(countrycode)

# worldbank income
worldbank = fread("~/Dropbox/nCoV/WorldBankIncome.csv")
worldbank[, iso2c := countrycode(Code, "wb", "iso2c")]

# get cities etc
data(world.cities) # NOTE THIS IS FROM 2006
wc = data.table(world.cities)
wc = wc[capital == 1]
wc[, cc := countrycode(country.etc, "country.name", "iso2c")]
prem_countries = names(matrices)[5:156];
prem_ccs = countrycode(prem_countries, "country.name", "iso2c")
wc = wc[cc %in% prem_ccs];
wc = wc[pop != 32187] # remove first San Jose
wc = wc[pop != 42372] # remove second Nicosia
wc[name == "Guatemala", name := "Guatemala (city)"]
wc[name == "Panama", name := "Panama (city)"]
wc[name == "Singapore", name := "Singapore (city)"]

populations = cm_populations
populations[, cc := countrycode(name, "country.name", "iso2c")]
populations[name == "Italy | Milan", country_code := 5555]

age_info = NULL

for (row in 1:nrow(wc)) {
    thiscc = wc[row]$cc;
    if (thiscc == "AD") { next; }
    if (thiscc == "MC") { next; }
    f = populations[cc == thiscc & country_code < 1000, f] * wc[row, pop] / (sum(populations[cc == thiscc, f + m]) * 1000);
    m = populations[cc == thiscc & country_code < 1000, m] * wc[row, pop] / (sum(populations[cc == thiscc, f + m]) * 1000);
    ages = populations[cc == thiscc & country_code < 1000, age];
    pop = data.table(country_code = 5555, name = wc[row]$name, age = ages, f = f, m = m, location_type = 5, cc = thiscc)
    age_info = rbind(age_info, pop[, .(name = unique(name), mean_age = 5 * weighted.mean(.I, f+m) - 2.5, median_age = 5 * weighted.median(.I, f+m) - 2.5)]);
    populations = rbind(populations, pop)
}

populations[country_code == 5555, unique(.N), by = name][V1 != 21]
populations = populations[name != "Andorra la Vella" & name != "Monaco-Ville"]
wc = wc[name != "Andorra la Vella" & name != "Monaco-Ville"]

wc = merge(wc, worldbank, by.x = "cc", by.y = "iso2c", all.x = T)

# RUN CAPITALS
argv = commandArgs();
argc = length(argv);
alt_lower = argv[argc - 2] == "shift"
nsamp = as.numeric(argv[argc - 1]);
which_fIa = as.numeric(argv[argc]);

pop_regions = wc[, name]
mat_regions = prem_countries[match(wc$cc, prem_ccs)]
loc_names = pop_regions
incomes = wc[, `Income group`]
if (alt_lower) {
    lower_income = incomes %in% c("Lower middle income", "Low income")
} else {
    lower_income = rep(F, length(incomes))
}

thousands = populations[name %in% pop_regions, .(popK = sum(f+m)), by = name]
thousands = merge(thousands, data.table(name = pop_regions, location = loc_names))
thousands[, name := NULL]

add_lines = function(data, name, fIa, alt_lower, k)
{
    filename = path(paste0("4-capitals-", name, ifelse(alt_lower, "-lowincome", ""), "-fix-u-fIa", which_fIa, ".csv"));
    if (k == 1) {
        fwrite(data, filename, append = FALSE, col.names = T)
    } else {
        fwrite(data, filename, append = T, col.names = F)
    }
}



cm_populations = populations

k = 1
for (fIa in which_fIa) {
    print(fIa);
    for (i in 1:length(pop_regions)) {
        cat(pop_regions[i]);
        for (j in 1:nsamp) {
            cat(".");
            p = cm_translate_parameters(parameters(pop_regions[i], mat_regions[i], "cov", -1, fIa, 1, 1, lower_income[i]));
            dynCO = cm_backend_simulate(p)[[1]][[1]];
            R0s = data.table(virus = "cov", location = pop_regions[i], fIa = fIa, run = j - 1, R0 = cm_calc_R0(p, 1));
            gc()
            
            ###
            setDT(dynCO);
            dynCO[, run := 1]
            dynCO = melt(dynCO, id.vars = c("run", "t", "population", "group"), variable.name = "compartment", value.name = "value");
            ###

            full_results = cbind(dynCO, run = j - 1, virus = "cov", location = loc_names[i], fIa = fIa, schools = "open")
            res2 = full_results[compartment == "cases"]
            res2[group <= 4, age := "0-19"]
            res2[group > 4 & group <= 12, age := "20-59"]
            res2[group > 12, age := "60+"]
            res2 = res2[, .(cases = sum(value)), by = .(t, location, age)]
            
            distrib = cased(dynCO, j - 1, "cov", loc_names[i], fIa, "open")
            
            results = cbind(shrink(dynCO), run = j - 1, virus = "cov", location = loc_names[i], fIa = fIa, schools = "open")
            
            mean_age = rbind(
                cbind(what = "cases", virus = "cov", location = loc_names[i], fIa = fIa, run = j - 1,
                      dynCO[compartment == "cases", .(mean_age = 5 * weighted.mean(group, value) + 2.5), by = t]),
                cbind(what = "subclinical", virus = "cov", location = loc_names[i], fIa = fIa, run = j - 1,
                      dynCO[compartment == "subclinical", .(mean_age = 5 * weighted.mean(group, value) + 2.5), by = t])
            )
            
            add_lines(R0s, "R0s", fIa, alt_lower, k);
            add_lines(res2, "cases_age", fIa, alt_lower, k);
            add_lines(distrib, "distrib", fIa, alt_lower, k);
            add_lines(results, "results", fIa, alt_lower, k);
            add_lines(mean_age, "mean_age", fIa, alt_lower, k);
            
            k = k + 1
        }
        cat("\n");
    }
}

gather = function(fIa, alt_lower)
{
    fn = function(name, fIa, alt_lower) path(paste0("4-capitals-", name, ifelse(alt_lower, "-lowincome", ""), "-fix-u-fIa", fIa, ".csv"));
    R0s = fread(fn("R0s", fIa, alt_lower));
    cases_age = fread(fn("cases_age", fIa, alt_lower));
    distrib = fread(fn("distrib", fIa, alt_lower));
    results = fread(fn("results", fIa, alt_lower));
    mean_age = fread(fn("mean_age", fIa, alt_lower));
    
    cm_save(list(epi = results, cases_age = cases_age, dist = distrib, age_info = age_info, wc = wc, mean_age = mean_age, R0 = R0s),
        path(paste0("4-capitals", ifelse(alt_lower, "-lowincome", ""), "-fix-u-fIa", fIa, ".qs")))
}

gather(which_fIa, alt_lower)
