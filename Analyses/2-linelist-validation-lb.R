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
cm_path = "~/Documents/ncov_age/covidm/";
cm_force_rebuild = F;
cm_version = 1;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# get fIa from command line
argv = commandArgs();
argc = length(argv);
assumed_fIa = as.numeric(argv[argc - 1]);
sfit = argv[argc];

# load Lombardy data
matrices = readRDS(path("../all_matrices.rds"));
onsets = fread(path("../lombardy-onsets.csv"))
confirmations = fread(path("../lombardy-confirmations.csv"))
onsets[, date := ymd("2020-02-10") + days_since_feb10]
confirmations[, date := ymd("2020-02-10") + days_since_feb10]
onsets[, onsets := round(onsets)]
confirmations[, confirmations := round(confirmations)]

# load posteriors from fitting
if (sfit == "ind") {
    fitted_symp_ind = qread(path(paste0("2-linelist_mn_both_fit_ind_fIa", assumed_fIa, "-rbzvih.qs")));
    fitted_symp_mean = unname(unlist(fitted_symp_ind[location == "Lombardia", lapply(.SD, mean), .SDcols = y_00:y_70]))
    fitted_symp_mean = rep(fitted_symp_mean, each = 2)
    fitted_susc_ind = qread(path(paste0("2-linelist_mn_both_fit_ind_fIa", assumed_fIa, "-rbzvih.qs")));
    fitted_susc_mean = unname(unlist(fitted_susc_ind[location == "Lombardia", lapply(.SD, mean), .SDcols = u_00:u_70]))
    fitted_susc_mean = rep(fitted_susc_mean, each = 2)
} else if (sfit == "con") {
    fitted_symp = qread(path(paste0("2-linelist_both_fit_fIa", assumed_fIa, "-rbzvih.qs")));
    fitted_symp_mean = unname(unlist(fitted_symp[, lapply(.SD, mean), .SDcols = y_00:y_70]))
    fitted_symp_mean = rep(fitted_symp_mean, each = 2)
    fitted_susc = qread(path(paste0("2-linelist_both_fit_fIa", assumed_fIa, "-rbzvih.qs")));
    fitted_susc_mean = unname(unlist(fitted_susc[, lapply(.SD, mean), .SDcols = u_00:u_70]))
    fitted_susc_mean = rep(fitted_susc_mean, each = 2)
} else {
    stop("Last argument must be ind or con.");
}

#
# LB VALIDATION - ALL IN ONE
#

# base parameters
lb_base_parameters = function(vary, fIa, report)
{
    # Preparing country-specific parameters
    mats = matrices[["Italy | Lombardia"]];
    n_age_groups = nrow(mats$home);
    size = cm_make_population("Italy | Lombardia", n_age_groups);
    
    # set global parameters
    p = list();
    p$time_step = 0.25;
    p$date0 = "2020-01-01";
    p$time0 = "2020-01-01";
    p$time1 = "2020-03-15";
    p$report_every = 4;
    p$travel = matrix(1, nrow = 1, ncol = 1);
    p$fast_multinomial = F;
    p$deterministic = T;
    
    # set population parameters
    pop = list();
    pop$type = "SEI3R";
    
    pop$size = size;
    pop$matrices = mats;
    pop$contact = c(1, 1, 1, 1);
    pop$contact_mult = c(1, 1, 1, 1);
    pop$contact_lowerto = c(100, 100, 100, 100);

    pop$dE  = cm_delay_gamma(3.0, 4.0, t_max = 60, t_step = 0.25)$p # Derived from He et al (44% infectiousness presymptomatic) https://www.nature.com/articles/s41591-020-0869-5#Sec9
    pop$dIp = cm_delay_gamma(2.1, 4.0, t_max = 60, t_step = 0.25)$p # Derived from Lauer et al (5.1 day incubation period) https://www.ncbi.nlm.nih.gov/pubmed/32150748
    pop$dIs = cm_delay_gamma(2.9, 4.0, t_max = 60, t_step = 0.25)$p # 5 days total: 5.5 day serial interval
    pop$dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p # Assumed same infectious period as clinical cases
    pop$dH  = c(1, 0); # hospitalization ignored
    pop$dC  = cm_delay_gamma(7.0, 4.0, t_max = 60, t_step = 0.25)$p # estimated below
    
    pop$fIp = rep(1, n_age_groups);
    pop$fIs = rep(1, n_age_groups);
    pop$fIa = rep(fIa, n_age_groups);
    pop$rho = rep(report, n_age_groups);
    pop$tau = rep(1, n_age_groups);
    pop$seed_times = c()
    pop$dist_seed_ages = rep(1, n_age_groups);
    pop$schedule = list();
    
    p$pop[[1]] = pop;

    return (p)
}

# Evaluate log-likelihood of model fit
likelihood_lb = function(parameters, dynamics, data, theta)
{
    theta = as.list(theta)
    
    ### Evaluate case onsets, just up to March 1st
    onsets2 = onsets[date <= ymd("2020-03-01")]
    onsets2[, t := as.numeric(date - ymd(parameters$date0))];
    y = merge(dynamics[compartment == "cases", .(model_case = sum(value)), by = t],
              onsets2, by = "t");
    ll_onsets = sum(dnbinom(y$onsets, size = 10, mu = pmax(0.1, 0.5 * theta$onset_known * y$model_case), log = T));
    
    ### Evaluate case confirmations
    conf2 = confirmations[date <= ymd("2020-03-01")]
    conf2[, t := as.numeric(date - ymd(parameters$date0))];
    y = merge(dynamics[compartment == "cases_reported", .(model_rep = sum(value)), by = t],
              conf2, by = "t");
    ll_confirmations = sum(dnbinom(y$confirmations, size = 10, mu = pmax(0.1, y$model_rep), log = T));

    return (ll_onsets + ll_confirmations)
}


### FITTING ###
pf_symp = function(p, x)
{
    x = as.list(x);
    n_age_groups = nrow(p$pop[[1]]$matrices$home);

    p$pop[[1]]$u = x$susc * c(fitted_susc_mean, rep(tail(fitted_susc_mean, 1), n_age_groups - length(fitted_susc_mean)));
    p$pop[[1]]$y = c(fitted_symp_mean, rep(tail(fitted_symp_mean, 1), n_age_groups - length(fitted_symp_mean)));
    p$pop[[1]]$dC  = cm_delay_gamma(x$conf_mean, x$conf_shape, t_max = 60, t_step = 0.25)$p

    p$pop[[1]]$dist_seed_ages = cm_age_coefficients(20, 50, 5 * (0:n_age_groups));
    p$pop[[1]]$seed_times = round(x$seed_t0):round(x$seed_t0 + x$seed_d);
    p$pop[[1]]$schedule = list(
        list(t = round(x$lockdown_t), contact = c(x$qL, x$qL, x$qL, x$qL))
    );

    return (p);
}


priors = list(
    susc = "N 0.1 0.025 T 0 0.25",
    conf_mean = "N 7 10 T 0 21",
    conf_shape = "N 7 10 T 0 21",
    onset_known = "B 2 2",
    seed_t0 = "N 30 30 T 0 60", 
    seed_d = "B 2 2 S 0 7 T 0 7",
    lockdown_t = "U 36 76",
    qL = "B 2 2"
);

pf = pf_symp;

fit = cm_fit(
    base_parameters = lb_base_parameters("symp", assumed_fIa, 0.2),
    priors = priors,
    parameters_func = pf,
    likelihood_func = likelihood_lb,
    data = 0,
    mcmc_burn_in = 3000, mcmc_samples = 10000, mcmc_init_opt = T, opt_local = F, opt_global_maxeval = 5000
);
cm_save(fit, path(paste0("2-linelist-validation-both-lb-fIa-", assumed_fIa, "-", sfit, ".qs")));


