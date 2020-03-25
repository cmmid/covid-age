# Load packages
library(Rcpp) 
library(RcppGSL)
library(HDInterval)
library(ggplot2)
library(data.table)
library(socialmixr)
library(shiny)
library(lubridate)
library(readxl)
library(cowplot)
library(stringr)
library(Hmisc)
library(extraDistr)
library(nloptr)
library(viridis)
library(magick)
library(qs)

hdimean = function(x, credMass = 0.95)
{
    h = hdi(x, credMass);
    m = mean(x);
    return (list(mean = m, lower = h[[1]], upper = h[[2]]))
}

htitle = function(text, size = 6, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5) {
    ggdraw() + draw_label(text, x = x, y = y, size = size, hjust = hjust, vjust = vjust)
}

# Load data
path = function(x) paste0("~/Dropbox/nCoV/", x);
matrices = readRDS(path("all_matrices.rds"));
populations = readRDS(path("wpp2019_pop2020.rds"));
age_dist = fread(path("Analyses/age.csv"))

# Fitting interface
sourceCpp(path("Corona/RMCMC/RMCMC.cpp"), verbose = T);
RMCMC = function(...) { return (data.table(MCMC(...))) }
ROptimize = function(...) { return (as.data.table(Optimize(...))) }

make_data = function(loc)
{
    if (loc == "Singapore") {
        pop.location = "Singapore";
        mat.location = "Singapore";
        mat.h = 1; mat.w = 1; mat.s = 1; mat.o = 1; # Considering Singapore's measures comprehensive across the board.
        weight = 1;
    } else if (loc == "South Korea") {
        pop.location = "Republic of Korea";
        mat.location = "Republic of Korea";
        mat.h = 1; mat.w = 1; mat.s = 1; mat.o = 1; # Generally open, then closed from late Feb. # TODO CHECK THIS. WERE SCHOOLS MOSTLY CLOSED?
        weight = 1;
    } else if (loc == "Japan") {
        pop.location = "Japan";
        mat.location = "Japan";
        mat.h = 1; mat.w = 1; mat.s = 1; mat.o = 1; # Schools open until March 2nd.
        weight = 1;
    } else if (loc == "Wuhan_CCDC") {
        pop.location = "China | Wuhan";
        mat.location = "China | Wuhan";
        mat.h = 1; mat.w = 1; mat.s = 0; mat.o = 1; # schools closed for most of outbreak
        weight = 1;
    } else if (loc == "Hubei_CCDC") {
        pop.location = "China | Hubei";
        mat.location = "China | Hubei";
        mat.h = 1; mat.w = 1; mat.s = 0; mat.o = 1; # schools closed for most of outbreak
        weight = 1;
    } else if (loc == "ChinaNonHubei_CCDC") {
        pop.location = "China";
        mat.location = "China | China";
        mat.h = 1; mat.w = 1; mat.s = 0; mat.o = 1; # schools closed for most of outbreak
        weight = 1;
    } else if (loc == "Ontario") {
        pop.location = "Canada";
        mat.location = "Canada";
        mat.h = 1; mat.w = 1; mat.s = 1; mat.o = 1; # Schools open until March 23rd
        weight = 1;
    } else if (loc %in% c("Anhui", "Guangdong", "Guangxi", "Hubei", "Hunan", 
        "Jiangsu", "Jiangxi", "Jilin", "Shaanxi", "Shandong", "Sichuan", "Tianjin", "Zhejiang")) { # province of China
        pop.location = paste("China |", loc);
        mat.location = paste("China |", loc);
        mat.h = 1; mat.w = 1; mat.s = 0; mat.o = 1; # schools closed for most of outbreak
        weight = 1/13;
    } else { # region of Italy
        pop.location = paste("Italy |", loc);
        mat.location = paste("Italy |", loc);
        mat.h = 1; mat.w = 1; mat.s = 1; mat.o = 1; # schools open for most of outbreak
        weight = 1/12;
    }

    # MATRIX
    m = matrices[[mat.location]];
    mat = mat.h * m$home + mat.w * m$work + mat.s * m$school + mat.o * m$other;

    n_age_groups = nrow(mat);

    # POPULATION
    pop = populations[name == pop.location, f + m];
    if (length(pop) > n_age_groups) {
        pop = c(pop[1:(n_age_groups - 1)], sum(pop[n_age_groups:length(pop)]));
    } else if (length(pop) < n_age_groups) {
        stop("Length of population distribution is less than size of contact matrix.");
    }
    
    # DATA
    counts = age_dist[location == loc, .(age_lower, age_upper, n)]
    
    return (list(mat = mat, pop = pop, counts = counts, weight = weight))
}

make_priors = function(np, prior, prefix = "f_")
{
    p = as.list(rep(prior, np))
    names(p) = sprintf("%s%02d", prefix, (0:(np-1))*10);
    return (p)
}

# Details of fitting

# expected_case_dist
#  contact: n * n contact matrix (number of group-j people contacted by a group-i person per day)
#  pop: n-vector, population distribution
#  d_Ip, Is, Ia: duration of presymptomatic, symptomatic, asymptomatic states
#  f_Ip, Is, Ia: relative infectiousness of presymptomatic, symptomatic, asymptomatic states
#  susc: n-vector, susceptibility of group-i individual
#  symp: n-vector, symptomatic probability of group-i individual
# returns n-vector, expected distribution of symptomatic cases across age groups; sums to 1
expected_case_dist = function(contact, pop, d_Ip, d_Is, d_Ia, f_Ip, f_Is, f_Ia, susc, symp)
{
    ngm = susc * t(t(contact) * (
        symp * (f_Ip * d_Ip + f_Is * d_Is) + 
        (1 - symp) * f_Ia * d_Ia)
    )
    dist = symp * pop * abs(eigen(ngm)$vectors[,1])
    dist / sum(dist)
}

# likelihood of model fit
#  x: parameters for model fit
#  d: list containing
#  d$vary: "susc", "symp", or "both"
#  d$symp: value of symp when not being fit
#  d$mode: NOT USED
#  d$data: data.frame with columns age_min, age_max, n
#  d$contact: contact matrix
#  d$pop: population vector
#  d$d_Ip
#  d$d_Is
#  d$d_Ia
#  d$f_Ip
#  d$f_Is
#  d$f_Ia
likelihood = function(x, D)
{
    ll_data = 0;
    if (D$method != "multinomial") {
        size = x[length(x)];
        x = x[1:(length(x) - 1)];
    }
    
    for (i in 1:length(D$data)) {
        
        d = D$data[[i]];
        
        # sanity check on d parameter
        n_age_groups = nrow(d$mat);
        if (ncol(d$mat) != n_age_groups | length(d$pop) != n_age_groups) {
            stop("likelihood: contact must be an NxN matrix, and pop an N-vector.");
        }
        
        # construct susc and symp
        if (D$vary == "susc") {
            susc = c(rep(x, each = 2), rep(x[length(x)], n_age_groups - length(x)*2));
            path = x;
            symp = D$symp;
        } else if (D$vary == "symp") {
            susc = rep(1, n_age_groups);
            symp = c(rep(x, each = 2), rep(x[length(x)], n_age_groups - length(x)*2));
            path = x;
        } else {
            stop("likelihood: vary must be 'susc' or 'symp'.");
        }
        
        # expected case distribution given model parameters
        dist = expected_case_dist(d$mat, d$pop, 
            D$d_Ip, D$d_Is, D$d_Ia, 
            D$f_Ip, D$f_Is, D$f_Ia,
            susc, symp);
    
        # amalgamate into 10-year bands
        dist = dist[seq(1, length(dist), by = 2)] + dist[seq(2, length(dist), by = 2)];
        counts = d$counts$n;
        if (length(counts) > length(dist)) {
            counts = c(counts[1:(length(dist)-1)], sum(counts[length(dist):length(counts)]));
        }
    
        # likelihood of data given model
        if (D$method == "multinomial") {
            ll_data = ll_data + sum(dmultinom(counts, sum(counts), dist, log = T)) * d$weight;
        } else {
            ll_data = ll_data + sum(ddirmnom(counts, sum(counts), size * dist, log = T)) * d$weight;
        }
    }

    ll_gauss = sum(dnorm(diff(path), 0, D$g_sd, log = T));
    
    return (ll_data + ll_gauss)
}

# Do fitting
do_fit = function(age_dist, varying, locations, method, fIa)
{
    # Get data for fitting
    fit = NULL;
    d = list(
        vary = varying,
        symp = 0.5,
        data = lapply(locations, make_data),
        d_Ip = 2.4,#
        d_Is = 3.2,#
        d_Ia = 7,#
        f_Ip = 1,
        f_Is = 1,
        f_Ia = fIa,
        g_sd = 0.25,
        method = method
    );

    # Adjust weights to multiply to 1    
    wt = 0;
    for (i in 1:length(d$data)) {
        wt = wt + log(d$data[[i]]$weight);
    }
    for (i in 1:length(d$data)) {
        d$data[[i]]$weight = d$data[[i]]$weight / exp(wt / length(d$data));
    }

    # Set up priors
    if (method == "multinomial") {
        priors = make_priors(8, "B 2 2");
    } else {
        priors = c(make_priors(8, "B 2 2"), size = "N 0 100 T 0 1000");
    }
    
    fit = RMCMC(likelihood, d, 2500, 1000, 18, priors);

    return (fit)
}

all_locations = age_dist[, unique(location)];

argv = commandArgs();
argc = length(argv);
fIa = as.numeric(argv[argc]);

f_ind = NULL;
for (loc in all_locations) {
    cat(paste0(loc, "...\n"));
    f = do_fit(age_dist, "symp", loc, "multinomial", fIa);
    f_ind = rbind(f_ind, cbind(f, location = loc))
}
qsave(f_ind, path(paste0("Analyses/2-linelist_symp_fit_ind_fIa", fIa, ".qs")), "balanced");

f_all = do_fit(age_dist, "symp", all_locations, "dirichlet-multinomial", fIa)
qsave(f_all, path(paste0("Analyses/2-linelist_symp_fit_fIa", fIa, ".qs")), "balanced");

# make expected case distributions

get_expected = function(x, D)
{
    n_age_groups = nrow(D$data[[1]]$mat);
    
    # construct susc and symp
    if (D$vary == "susc") {
        susc = c(rep(x, each = 2), rep(x[length(x)], n_age_groups - length(x)*2));
        symp = D$symp;
    } else if (D$vary == "symp") {
        susc = rep(1, n_age_groups);
        symp = c(rep(x, each = 2), rep(x[length(x)], n_age_groups - length(x)*2));
    } else {
        stop("likelihood: vary must be 'susc' or 'symp'.");
    }
        
    # expected case distribution given model parameters
    dist = expected_case_dist(D$data[[1]]$mat, D$data[[1]]$pop, 
        D$d_Ip, D$d_Is, D$d_Ia, 
        D$f_Ip, D$f_Is, D$f_Ia,
        susc, symp);
    
    # amalgamate into 10-year bands
    dist = dist[seq(1, length(dist), by = 2)] + dist[seq(2, length(dist), by = 2)];
    if (length(dist) > length(x)) {
        dist = c(dist[1:(length(x)-1)], sum(dist[length(x):length(dist)]));
    }
    
    return (dist)
}

results = NULL;
for (loc in all_locations)
{
    d = list(
        vary = "symp",
        symp = 0.5,
        data = lapply(loc, make_data),
        d_Ip = 2.4,#
        d_Is = 3.2,#
        d_Ia = 7,#
        f_Ip = 1,
        f_Is = 1,
        f_Ia = fIa,
        g_sd = 0.15
    );

    draws = NULL
    data = f_ind[location == loc];
    data = data[sample.int(.N, 100)];
    for (i in 1:nrow(data))
    {
        cd = get_expected(unlist(unname(data[i, 5:12])), d);
        names(cd) = names(data[i, 5:12]);
        draws = rbind(draws, as.data.table(as.list(cd)));
    }
    
    results = rbind(results,
        cbind(draws, location = loc));
};

qsave(results, path(paste0("Analyses/2-linelist_symp_fit_ind_expected_fIa", fIa, ".qs")), "balanced")
