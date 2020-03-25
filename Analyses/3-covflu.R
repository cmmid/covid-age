library(ggplot2)
library(data.table)
library(lubridate)
library(stringr)
library(readxl)

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_builddir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

flu_scenario = data.table(read_excel(path("../flu_scenario_v3.xlsx")));
flu_scenario[, Susceptibility := OldSusc]
flu_scenario[, Severity := scaled_newSev]
covid_scenario = qread(path("2-linelist_symp_fit_fIa0.5.qs"));
variability = mean(unlist(covid_scenario[, lapply(.SD, function(x) sd(x)/mean(x)), .SDcols = f_00:f_70]));
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
parameters = function(pop_region, mat_region, virus, u, R0, fIa, report, qS, dyn_trigger = NULL, inc_trigger = NULL)
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
    p$time1 = "2021-01-01";
    p$report_every = 4;
    p$travel = matrix(1.0, nrow = 1, ncol = 1);
    p$fast_multinomial = F;
    p$deterministic = T;
    
    # set population parameters
    pop = list();
    pop$type = "SEI3R";
    
    pop$dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p # Derived from Backer et al Eurosurveillance
    pop$dIp = cm_delay_gamma(2.4, 4.0, t_max = 60, t_step = 0.25)$p # Derived from Backer et al Eurosurveillance
    pop$dIa = cm_delay_gamma(7.0, 4.0, t_max = 60, t_step = 0.25)$p # Assumed 7 days subclinical shedding
    pop$dIs = cm_delay_gamma(3.2, 3.7, t_max = 60, t_step = 0.25)$p # Zhang et al 2020
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
        pop$u = 0.1 * c(flu_scenario$Susceptibility[1:(n_age_groups - 1)], mean(flu_scenario$Susceptibility[n_age_groups:nrow(flu_scenario)]));
        pop$y = c(flu_scenario$Severity[1:(n_age_groups - 1)], mean(flu_scenario$Severity[n_age_groups:nrow(flu_scenario)]));
        pop$y = pmax(0, pmin(1, pop$y * rnorm(length(pop$y), 1, variability)));
    } else {
        pop$u = rep(0.1, n_age_groups);
        covy = unname(unlist(covid_scenario[sample.int(nrow(covid_scenario), 1), f_00:f_70]));
        covy = unname(rep(covy, each = 2))
        if (length(covy) < n_age_groups) {
            covy = c(covy, rep(tail(covy, 1), n_age_groups - length(covy)));
        } else if (n_age_groups < length(covy)) {
            covy = c(covy[1:(n_age_groups - 1)], mean(covy[n_age_groups:length(covy)]));
        }
        pop$y = covy;
    }

    p$pop[[1]] = pop;
    
    if ((R0 >= 0 & u >= 0) | (R0 < 0 & u < 0)) {
        stop("Specify exactly one of R0 and u, setting the other to a negative number.");
    }
    
    if (R0 >= 0) {
        R0_actual = cm_calc_R0(p, 1);
        factor = R0 / R0_actual;
        pop$u = factor * pop$u;
    } else if (u >= 0) {
        if (virus == "flu") {
            pop$u = u * c(flu_scenario$Susceptibility[1:(n_age_groups - 1)], mean(flu_scenario$Susceptibility[n_age_groups:nrow(flu_scenario)]));
        } else {
            pop$u = rep(u, n_age_groups);
        }
    }

    if (!is.null(dyn_trigger)) {
        setDT(dyn_trigger)
        trig = dyn_trigger[compartment == "cases", .(cases = sum(value)), by = t];
        trig_time = head(trig[cases > inc_trigger, t], 1);
        pop$schedule = list(
            list(t = trig_time, contact = c(1, 1, qS, 1)),
            list(t = trig_time + 90, contact = c(1, 1, 1, 1))
        )
    }

    p$pop[[1]] = pop;
    
    return (cm_translate_parameters(p))
}

shrink = function(dyn)
{
    setDT(dyn)
    dyn[compartment == "cases", .(cases = sum(value)), by = t]
}

cased = function(dyn, run, virus, location, fIa, schools)
{
    setDT(dyn)
    cbind(cm_case_distribution(dyn, "2020-01-01", "2020-01-01", "2021-01-01", seq(0, 80, by = 10)),
          run = run, virus = virus, location = location, fIa = fIa, schools = schools)
}

pop_regions = c("UK | Birmingham", "Zimbabwe | Bulawayo", "Italy | Milan")
mat_regions = c("UK | Birmingham", "Zimbabwe | Bulawayo", "Italy | Milan")
loc_names = c("Birmingham", "Bulawayo", "Milan");

thousands = cm_populations[name %in% pop_regions, .(popK = sum(f+m)), by = name]
thousands = merge(thousands, data.table(name = pop_regions, location = loc_names))
thousands[, name := NULL]

argv = commandArgs();
argc = length(argv);
a_u = as.numeric(argv[argc - 1]);
a_R0 = as.numeric(argv[argc]);

print(a_u)
print(a_R0)

results = NULL;
distrib = NULL;
R0s = NULL;
nsamp = 100;
for (fIa in c(0, 0.25, 0.5, 0.75)) {
    print(fIa);
    covid_scenario = qread(path(paste0("2-linelist_symp_fit_fIa", fIa, ".qs")));
    variability = mean(unlist(covid_scenario[, lapply(.SD, function(x) sd(x)/mean(x)), .SDcols = f_00:f_70]));
    
    if (a_R0 < 0) {
        if (fIa == 0) {
            u_flu = 0.985;
            u_cov = 0.077;
            R0_flu = -1;
            R0_cov = -1;
        } else if (fIa == 0.25) {
            u_flu = 0.148;
            u_cov = 0.053;
            R0_flu = -1;
            R0_cov = -1;
        } else if (fIa == 0.5) {
            u_flu = 0.076;
            u_cov = 0.039;
            R0_flu = -1;
            R0_cov = -1;
        } else {
            u_flu = 0.051;
            u_cov = 0.03;
            R0_flu = -1;
            R0_cov = -1;
        }
    } else {
        u_flu = -1;
        u_cov = -1;
        R0_flu = 2.4;
        R0_cov = 2.4;
    }

    for (i in 1:length(pop_regions)) {
        cat(pop_regions[i]);
        for (j in 1:nsamp) {
            cat(".");
            p = cm_translate_parameters(parameters(pop_regions[i], mat_regions[i], "flu", u_flu, R0_flu, fIa, 1, 1));
            R0s = rbind(R0s, data.table(virus = "flu", a_u = a_u, a_R0 = a_R0, fIa = fIa, location = loc_names[i], R0 = cm_calc_R0(p, 1)));
            dynFO = cm_backend_simulate(p);
            setDT(dynFO);
            dynFO = melt(dynFO, id.vars = c("run", "t", "population", "group"), variable.name = "compartment", value.name = "value");

            p = cm_translate_parameters(parameters(pop_regions[i], mat_regions[i], "flu", u_flu, R0_flu, fIa, 1, 0, dynFO, thousands[location == loc_names[i], popK] * 0.01))
            dynFC = cm_backend_simulate(p);
            setDT(dynFC);
            dynFC = melt(dynFC, id.vars = c("run", "t", "population", "group"), variable.name = "compartment", value.name = "value");

            p = cm_translate_parameters(parameters(pop_regions[i], mat_regions[i], "cov", u_cov, R0_cov, fIa, 1, 1));
            R0s = rbind(R0s, data.table(virus = "cov", a_u = a_u, a_R0 = a_R0, fIa = fIa, location = loc_names[i], R0 = cm_calc_R0(p, 1)));
            dynCO = cm_backend_simulate(p);
            setDT(dynCO);
            dynCO = melt(dynCO, id.vars = c("run", "t", "population", "group"), variable.name = "compartment", value.name = "value");

            p = cm_translate_parameters(parameters(pop_regions[i], mat_regions[i], "cov", u_cov, R0_cov, fIa, 1, 0, dynCO, thousands[location == loc_names[i], popK] * 0.01));
            dynCC = cm_backend_simulate(p);
            setDT(dynCC);
            dynCC = melt(dynCC, id.vars = c("run", "t", "population", "group"), variable.name = "compartment", value.name = "value");
            gc()
            
            distrib = rbind(distrib,
                cased(dynFO, j - 1, "flu", loc_names[i], fIa, "open"),
                cased(dynFC, j - 1, "flu", loc_names[i], fIa, "closed"),
                cased(dynCO, j - 1, "cov", loc_names[i], fIa, "open"),
                cased(dynCC, j - 1, "cov", loc_names[i], fIa, "closed")
            );
            
            results = rbind(results,
                cbind(shrink(dynFO), run = j - 1, virus = "flu", location = loc_names[i], fIa = fIa, schools = "open"),
                cbind(shrink(dynFC), run = j - 1, virus = "flu", location = loc_names[i], fIa = fIa, schools = "closed"),
                cbind(shrink(dynCO), run = j - 1, virus = "cov", location = loc_names[i], fIa = fIa, schools = "open"),
                cbind(shrink(dynCC), run = j - 1, virus = "cov", location = loc_names[i], fIa = fIa, schools = "closed")
            );
        }
        cat("\n");
    }
}

results = merge(results, thousands, by = "location")

results[, schools := factor(schools, levels = c("open", "closed"))]
results[, location := factor(location, levels = c("Milan", "Birmingham", "Bulawayo"))]

if (a_u < 0) {
    cm_save(list(epi = results, dist = distrib, R0 = R0s), path("3-covflu-fix-R0.qs"))
} else if (a_R0 < 0) {
    cm_save(list(epi = results, dist = distrib, R0 = R0s), path("3-covflu-fix-u.qs"))
}
