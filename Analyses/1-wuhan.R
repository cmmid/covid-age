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

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Documents/ncov_age/covidm/";
cm_version = 1;
cm_force_rebuild = F;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

age_dist_wuhan = fread("~/Dropbox/nCoV/wuhan_dist_70.csv");
age_dist_wuhan[, age := factor(age, levels = unique(age))];
ciw = cbind(age_dist_wuhan[, rdirichlet(1000, n + 1) * sum(n), by = date][, .(date, x = V1)], age_dist_wuhan[, rep(age, each = 1000), by = date][, .(age = V1)]);
ciw = ciw[, cm_mean_hdi(x), keyby = .(date, age)]
age_dist_wuhan = merge(age_dist_wuhan, ciw, by = c("date", "age"))
age_dist_wuhan = age_dist_wuhan[, .(age = age, mean = n/sum(n), lower = lower/sum(n), upper = upper / sum(n)), by = date]

cases_wuhan = fread("~/Dropbox/nCoV/44000_incidence.csv")
cases_wuhan[, date := ymd(date)]
cases_wuhan[, incidence := round(incidence)]
cases_wuhan1 = cases_wuhan[date <= ymd("2020-01-24")]
cases_wuhan2 = data.table(date0 = ymd("2020-01-25"), date1 = ymd("2020-02-01"))
cases_wuhan2$incidence = cases_wuhan[date >= cases_wuhan2$date0 & date <= cases_wuhan2$date1, sum(incidence)]

# FITTING

# base parameters
base_parameters = function(pop_region, mat_region, susc, symp, fIa, report)
{
    # p = cm_parameters_SEI3R(pop_region, mat_region, date_start = "2019-11-01", date_end = "2020-02-14");
    # 
    # p$pop[[1]]$dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p # Derived from Backer et al Eurosurveillance
    # p$pop[[1]]$dIp = cm_delay_gamma(2.4, 4.0, t_max = 60, t_step = 0.25)$p # Derived from Backer et al Eurosurveillance
    # p$pop[[1]]$dIa = cm_delay_gamma(7.0, 4.0, t_max = 60, t_step = 0.25)$p # Assumed 7 days subclinical shedding
    # p$pop[[1]]$dIs = cm_delay_gamma(3.2, 3.7, t_max = 60, t_step = 0.25)$p # Zhang et al 2020
    # p$pop[[1]]$dH  = c(1, 0); # hospitalization ignored
    # p$pop[[1]]$dC  = c(1, 0); # cases are reported immediately -- because fitting is retrospective to case onset.
    # 
    # n_age_groups = length(p$pop[[1]]$size);
    # 
    # p$pop[[1]]$u = rep(susc, n_age_groups);
    # p$pop[[1]]$y = rep(symp, n_age_groups);
    # p$pop[[1]]$fIa = rep(fIa, n_age_groups)
    # p$pop[[1]]$rho = rep(report, n_age_groups)
    # 
    # return (p);
    

    # Preparing country-specific parameters
    mats = cm_matrices[[mat_region]];
    n_age_groups = nrow(mats$home);
    pop_size = cm_make_population(pop_region, n_age_groups);
    
    # set global parameters
    p = list();
    p$time_step = 0.25;
    p$date0 = "2019-11-01";
    p$time0 = "2019-11-01";
    p$time1 = "2020-02-14";
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
    pop$u = rep(susc, n_age_groups);
    pop$y = rep(symp, n_age_groups);
    pop$fIp = rep(1, n_age_groups);
    pop$fIs = rep(1, n_age_groups);
    pop$fIa = rep(fIa, n_age_groups)
    pop$rho = rep(report, n_age_groups);
    pop$tau = rep(1, n_age_groups);
    pop$seed_times = c();
    pop$dist_seed_ages = rep(1, n_age_groups);
    pop$schedule = list()

    p$pop[[1]] = pop;

    return (p)
}

# Evaluate log-likelihood of model fit
likelihood_wuhan = function(parameters, dynamics, data, theta)
{
    sz = 200
    ssum = function(x) { ifelse(sum(x) == 0, 1, sum(x)) }
    
    ### Evaluate incidence likelihood, pt.1
    cases_wuhan1[, t := as.numeric(date - ymd(parameters$date0))]
    incidence1 = merge(dynamics[compartment == "cases_reported", .(model_case = sum(value)), by = t],
              cases_wuhan1, by = "t");
    ll_incidence1 = sum(dnbinom(incidence1$incidence, size = sz, mu = pmax(0.1, incidence1$model_case), log = T));
    #ll_incidence1 = sum(dpois(incidence1$incidence, lambda = pmax(0.1, incidence1$model_case), log = T));

    ### Evaluate incidence likelihood, pt.2
    cases_wuhan2[, t0 := as.numeric(date0 - ymd(parameters$date0))]
    cases_wuhan2[, t1 := as.numeric(date1 - ymd(parameters$date0))]
    incidence2 = dynamics[compartment == "cases_reported" & t >= cases_wuhan2$t0 & t <= cases_wuhan2$t1, sum(value)];
    ll_incidence2 = sum(dnbinom(cases_wuhan2$incidence, size = sz, mu = pmax(0.1, incidence2), log = T));
    #ll_incidence2 = sum(dpois(cases_wuhan2$incidence, lambda = pmax(0.1, incidence2), log = T));
    
    ### Evaluate age distribution likelihood, pt. 1
    age_dist = dynamics[, cm_case_distribution(.SD, parameters$date0, "2019-12-08", c("2019-12-31", "2020-01-11", "2020-01-22"), c(0, 15, 45, 65, 100), "cases_reported")];
    age_dist[, fcases := pmax(1e-6, fcases)];
    age_dist[, measured := c(0, 12, 24, 11,
                            0, 39, 106, 103,
                            0, 33, 49, 48)];
    ll_agedist1 = age_dist[, ddirmnom(measured, sum(measured), sz * pmax(1e-4, fcases)/ssum(fcases), log = T), by = date][, sum(V1)];
    #ll_agedist1 = age_dist[, dmultinom(measured, sum(measured), fcases/sum(fcases), log = T), by = date][, sum(V1)];

    ### Evaluate age distribution likelihood, pt. 2
    age_dist = dynamics[, cm_case_distribution(.SD, parameters$date0, "2019-12-08", "2020-02-11", c(seq(0, 70, by = 10), 100), "cases_reported")];
    age_dist[, fcases := pmax(1e-6, fcases)];
    age_dist[, measured := round(19558 * c(0.4, 0.4, 4.5, 13.1, 15.6, 22.0, 26.5, 17.6) / 100)]; # China CDC Weekly
    ll_agedist2 = age_dist[, ddirmnom(measured, sum(measured), sz * pmax(1e-4, fcases)/ssum(fcases), log = T)];
    #ll_agedist2 = age_dist[, dmultinom(measured, sum(measured), fcases/sum(fcases), log = T)];
    
    return (ll_incidence1 + ll_incidence2 + ll_agedist1 + ll_agedist2)
}


### FITTING ###
pf_same = function(p, x)
{
    x = as.list(x);
    n_age_groups = nrow(p$pop[[1]]$matrices$home);
    p$pop[[1]]$u = rep(x$susc, n_age_groups);
    p$pop[[1]]$dist_seed_ages = cm_age_coefficients(x$seed_age - x$seed_age_range, x$seed_age + x$seed_age_range, 5 * (0:n_age_groups));
    seed_start = floor(x$seed_start);
    p$pop[[1]]$seed_times = seed_start + 0:13; # seed for 14 days
    p$pop[[1]]$schedule = list(
        list(t = "2020-01-12", contact = c(x$qH, x$qH, 0, x$qH)),
        list(t = "2020-01-23", contact = c(x$qL, x$qL, 0, x$qL))
    );
    return (p);
}

pf_susc = function(p, x)
{
    x = as.list(x);
    n_age_groups = nrow(p$pop[[1]]$matrices$home);

    # definition of "young", "middle", and "old"
    young  = cm_interpolate_cos(seq(2.5, n_age_groups * 5 - 2.5, by = 5), x$age_y, 1, x$age_m, 0);
    old    = cm_interpolate_cos(seq(2.5, n_age_groups * 5 - 2.5, by = 5), x$age_m, 0, x$age_o, 1);
    middle = 1 - young - old;

    p$pop[[1]]$u = young * x$susc_y + middle * x$susc_m + old * x$susc_o;
    p$pop[[1]]$dist_seed_ages = cm_age_coefficients(x$seed_age - x$seed_age_range, x$seed_age + x$seed_age_range, 5 * (0:n_age_groups));
    seed_start = floor(x$seed_start);
    p$pop[[1]]$seed_times = seed_start + 0:13; # seed for 14 days
    p$pop[[1]]$schedule = list(
        list(t = "2020-01-12", contact = c(x$qH, x$qH, 0, x$qH)),
        list(t = "2020-01-23", contact = c(x$qL, x$qL, 0, x$qL))
    );
    return (p);
}

pf_symp = function(p, x)
{
    x = as.list(x);
    n_age_groups = nrow(p$pop[[1]]$matrices$home);

    # definition of "young", "middle", and "old"
    young  = cm_interpolate_cos(seq(2.5, n_age_groups * 5 - 2.5, by = 5), x$age_y, 1, x$age_m, 0);
    old    = cm_interpolate_cos(seq(2.5, n_age_groups * 5 - 2.5, by = 5), x$age_m, 0, x$age_o, 1);
    middle = 1 - young - old;

    p$pop[[1]]$u = rep(x$susc, n_age_groups);
    p$pop[[1]]$y = young * x$symp_y + middle * x$symp_m + old * x$symp_o;
    p$pop[[1]]$dist_seed_ages = cm_age_coefficients(x$seed_age - x$seed_age_range, x$seed_age + x$seed_age_range, 5 * (0:n_age_groups));
    seed_start = floor(x$seed_start);
    p$pop[[1]]$seed_times = seed_start + 0:13; # seed for 14 days
    p$pop[[1]]$schedule = list(
        list(t = "2020-01-12", contact = c(x$qH, x$qH, 0, x$qH)),
        list(t = "2020-01-23", contact = c(x$qL, x$qL, 0, x$qL))
    );
    return (p);
}

pf_both = function(p, x)
{
    x = as.list(x);
    n_age_groups = nrow(p$pop[[1]]$matrices$home);

    # definition of "young", "middle", and "old"
    young  = cm_interpolate_cos(seq(2.5, n_age_groups * 5 - 2.5, by = 5), x$age_y, 1, x$age_m, 0);
    old    = cm_interpolate_cos(seq(2.5, n_age_groups * 5 - 2.5, by = 5), x$age_m, 0, x$age_o, 1);
    middle = 1 - young - old;

    p$pop[[1]]$u = x$susc * (young * x$eff_y + middle * x$eff_m + old * x$eff_o);
    p$pop[[1]]$y =    0.5 * (young * x$eff_y + middle * x$eff_m + old * x$eff_o);
    p$pop[[1]]$dist_seed_ages = cm_age_coefficients(x$seed_age - x$seed_age_range, x$seed_age + x$seed_age_range, 5 * (0:n_age_groups));
    seed_start = floor(x$seed_start);
    p$pop[[1]]$seed_times = seed_start + 0:13; # seed for 14 days
    p$pop[[1]]$schedule = list(
        list(t = "2020-01-12", contact = c(x$qH, x$qH, 0, x$qH)),
        list(t = "2020-01-23", contact = c(x$qL, x$qL, 0, x$qL))
    );
    return (p);
}


run_fit = function(scenario, fIa, rho)
{
    #                         pop_region,      mat_region, susc, symp, fIa, rep
    bp = base_parameters("China | Wuhan", "China | Wuhan",  0.1,  0.5, fIa, rho);
    bp$deterministic = T;
    bp$fast_multinomial = F;
    
    if (scenario == "same") {
        priors = list(
            susc = "N 0.1 0.025 T 0 1",
            seed_start = "N 15 30 T 0 30",
            seed_age = "N 60 20 T 30 80", 
            seed_age_range = "B 2 2 S 0 10", 
            qH = "B 2 2 S 0 2", 
            qL = "B 2 2"
        );
        pf = pf_same;
    } else if (scenario == "susc") {
        priors = list(
            age_y = "N 15 15 T 0 30", 
            age_m = "N 45 15 T 30 60", 
            age_o = "N 75 15 T 60 90",
            susc_y = "N 0.1 0.025 T 0 0.25",
            susc_m = "N 0.1 0.025 T 0 0.25",
            susc_o = "N 0.1 0.025 T 0 0.25",
            seed_start = "N 15 30 T 0 30",
            seed_age = "N 60 20 T 30 80", 
            seed_age_range = "B 2 2 S 0 10", 
            qH = "B 2 2 S 0 2", 
            qL = "B 2 2"
        );
        pf = pf_susc;
    } else if (scenario == "symp") {
        priors = list(
            age_y = "N 15 15 T 0 30", 
            age_m = "N 45 15 T 30 60", 
            age_o = "N 75 15 T 60 90",
            susc = "N 0.1 0.025 T 0 0.25",
            symp_y = "N 0.5 0.1 T 0 0.5",
            symp_m = "N 0.5 0.1 T 0 1",
            symp_o = "N 0.5 0.1 T 0.5 1",
            seed_start = "N 15 30 T 0 30",
            seed_age = "N 60 20 T 30 80", 
            seed_age_range = "B 2 2 S 0 10", 
            qH = "B 2 2 S 0 2", 
            qL = "B 2 2"
        );
        pf = pf_symp;
    } else if (scenario == "both") {
        priors = list(
            age_y = "N 15 15 T 0 30", 
            age_m = "N 45 15 T 30 60", 
            age_o = "N 75 15 T 60 90",
            susc = "N 0.1 0.025 T 0 0.25",
            eff_y = "N 1 0.20 T 0 1",
            eff_m = "N 1 0.20 T 0.5 1.5",
            eff_o = "N 1 0.20 T 1 2",
            seed_start = "N 15 30 T 0 30",
            seed_age = "N 60 20 T 30 80", 
            seed_age_range = "B 2 2 S 0 10", 
            qH = "B 2 2 S 0 2", 
            qL = "B 2 2"
        );
        pf = pf_both;
    }
    
    fit = cm_fit(
        base_parameters = bp,
        priors = priors,
        parameters_func = pf,
        likelihood_func = likelihood_wuhan,
        data = 0,
        mcmc_burn_in = 2000, mcmc_samples = 10000, mcmc_init_opt = T, opt_local = F, opt_global_maxeval = 5000
    );

    cm_save(fit, path(paste0("1-wuhan-fit-", scenario, "-rho-", rho,  "-fIa-", fIa, ".qs")));
    
    return (fit)
}

# fit_same = run_fit("same", 0.5)
# fit_susc = run_fit("susc", 0.5)
# fit_symp3 = run_fit("symp", 0.5)
# 
# hist(fit_symp$posterior$lp, breaks = 200)
# 
# rows = sample.int(fit_symp$posterior[, .N], replace = T, prob = fit_symp$posterior[, exp(lp)])
# fit_symp2 = duplicate(fit_symp)
# fit_symp2$posterior = fit_symp2$posterior[rows]
# 
# cm_plot_posterior(fit_symp3)

argv = commandArgs();
argc = length(argv);
run_fit(argv[argc-2], as.numeric(argv[argc-1]), as.numeric(argv[argc]))
