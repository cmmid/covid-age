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

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# get fIa from command line
argv = commandArgs();
argc = length(argv);
assumed_fIa = as.numeric(argv[argc - 1]);
sfit = argv[argc];

# Chinese provinces
cprovinces = c(
"Anhui",
"Guangdong",
"Guangxi",
"Hubei",
"Hunan",
"Jiangsu",
"Jiangxi",
"Jilin",
"Shaanxi",
"Shandong",
"Sichuan",
"Tianjin",
"Zhejiang")


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

# load posteriors from fitting
if (sfit == "ind") {
    # incidence is from CCDC...
    fitted_symp_ind = qread(path(paste0("2-linelist_symp_fit_ind_fIa", assumed_fIa, ".qs")));
    fitted_symp_mean2 = unname(unlist(fitted_symp_ind[location %like% "CCDC", lapply(.SD, mean), .SDcols = f_00:f_70]))
    fitted_symp_mean = rep(fitted_symp_mean2, each = 2)
} else if (sfit == "con") {
    fitted_symp = qread(path(paste0("2-linelist_symp_fit_fIa", assumed_fIa, ".qs")));
    fitted_symp_mean = unname(unlist(fitted_symp[, lapply(.SD, mean), .SDcols = f_00:f_70]))
    fitted_symp_mean = rep(fitted_symp_mean, each = 2)
} else {
    stop("Last argument must be ind or con.");
}

# FITTING

# base parameters
beishang_base_parameters = function(vary, fIa, report)
{
    # Preparing country-specific parameters
    matsB = cm_matrices[["China | Beijing"]];
    matsS = cm_matrices[["China | Shanghai"]];
    n_age_groups = nrow(matsB$home);
    sizeB = cm_make_population("China | Beijing", n_age_groups);
    sizeS = cm_make_population("China | Shanghai", n_age_groups);
    
    # set global parameters
    p = list();
    p$time_step = 0.25;
    p$date0 = "2019-12-01";
    p$time0 = "2019-12-01";
    p$time1 = "2020-02-29"
    p$report_every = 4;
    p$travel = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2);
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
    
    pop$fIp = rep(1, n_age_groups);
    pop$fIs = rep(1, n_age_groups);
    pop$fIa = rep(fIa, n_age_groups);
    pop$rho = rep(report, n_age_groups);
    pop$tau = rep(1, n_age_groups);
    pop$seed_times = c()
    pop$dist_seed_ages = rep(1, n_age_groups);
    pop$schedule = list();
    
    # Beijing-specific
    popB = duplicate(pop);
    popB$size = sizeB;
    popB$matrices = matsB;
    popB$contact = c(1, 1, 1, 1);
    popB$contact_mult = c(1, 1, 1, 1);
    popB$contact_lowerto = c(100, 100, 100, 100);

    # Shanghai-specific
    popS = duplicate(pop);
    popS$size = sizeS;
    popS$matrices = matsS;
    popS$contact = c(1, 1, 1, 1);
    popS$contact_mult = c(1, 1, 1, 1);
    popS$contact_lowerto = c(100, 100, 100, 100);

    p$pop[[1]] = popB;
    p$pop[[2]] = popS;

    return (p)
}

# Evaluate log-likelihood of model fit
likelihood_beishang = function(parameters, dynamics, data, theta)
{
    resultsB = dynamics[population == 1]
    resultsS = dynamics[population == 2]

    ### Evaluate incidence likelihood, Beijing
    cases_B[, t := as.numeric(date - ymd(parameters$date0))]
    x = merge(resultsB[compartment == "cases_reported", .(model_case = sum(value)), by = t],
              cases_B, by = "t");
    ll_incidenceB = sum(dnbinom(x$incidence, size = 10, mu = pmax(0.1, x$model_case), log = T));

    ### Evaluate incidence likelihood, Shanghai
    cases_S[, t := as.numeric(date - ymd(parameters$date0))]
    x = merge(resultsS[compartment == "cases_reported", .(model_case = sum(value)), by = t],
              cases_S, by = "t");
    ll_incidenceS = sum(dnbinom(x$incidence, size = 10, mu = pmax(0.1, x$model_case), log = T));

    return (ll_incidenceB + ll_incidenceS)
}


### FITTING ###
pf_symp = function(p, x)
{
    x = as.list(x);
    n_age_groups = nrow(p$pop[[1]]$matrices$home);

    p$pop[[1]]$u = rep(x$susc, n_age_groups);
    p$pop[[1]]$y = c(fitted_symp_mean, rep(tail(fitted_symp_mean, 1), n_age_groups - length(fitted_symp_mean)));

    p$pop[[1]]$dist_seed_ages = cm_age_coefficients(20, 50, 5 * (0:n_age_groups));
    p$pop[[1]]$seed_times = round(x$B_seed_t0):round(x$B_seed_t0 + x$seed_d);
    p$pop[[1]]$schedule = list(
        list(t = "2020-01-12",        contact = c(x$qH, x$qH, 0, x$qH)),
        list(t = round(x$lockdown_t), contact = c(x$qL, x$qL, 0, x$qL))
    );

    p$pop[[2]]$u = rep(x$susc, n_age_groups);
    p$pop[[2]]$y = c(fitted_symp_mean, rep(tail(fitted_symp_mean, 1), n_age_groups - length(fitted_symp_mean)));
    
    p$pop[[2]]$dist_seed_ages = cm_age_coefficients(20, 50, 5 * (0:n_age_groups));
    p$pop[[2]]$seed_times = round(x$S_seed_t0):round(x$S_seed_t0 + x$seed_d);
    p$pop[[2]]$schedule = list(
        list(t = "2020-01-12",        contact = c(x$qH, x$qH, 0, x$qH)),
        list(t = round(x$lockdown_t), contact = c(x$qL, x$qL, 0, x$qL))
    );

    return (p);
}


priors = list(
    susc = "N 0.1 0.025 T 0 0.25",
    B_seed_t0 = "N 0 30 T 0 60", 
    S_seed_t0 = "N 0 30 T 0 60", 
    seed_d = "B 2 2 S 0 7 T 0 7",
    lockdown_t = "U 43 100",
    qH = "B 2 2 S 0 2",
    qL = "B 2 2"
);

pf = pf_symp;

fit = cm_fit(
    base_parameters = beishang_base_parameters("symp", assumed_fIa, 0.2),
    priors = priors,
    parameters_func = pf,
    likelihood_func = likelihood_beishang,
    data = 0,
    mcmc_burn_in = 3000, mcmc_samples = 10000, mcmc_init_opt = T, opt_local = F, opt_global_maxeval = 5000
);
cm_save(fit, path(paste0("2-linelist-validation-bs-fIa-", assumed_fIa, "-", sfit, ".qs")));
