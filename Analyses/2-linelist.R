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

# Load data
path = function(x) paste0("~/Dropbox/nCoV/", x);
matrices = readRDS(path("all_matrices.rds"));
populations = readRDS(path("wpp2019_pop2020.rds"));
age_dist = fread(path("Analyses/age.csv"))

# Fitting interface
cm_path = "~/Documents/ncov_age/covidm/";
cm_force_rebuild = F;
if (Sys.info()["nodename"] %like% "lshtm") {
    cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(paste0(cm_path, "R/covidm.R"))
RMCMC = function(...) { data.table(cm_backend_mcmc(...)); }
ROptimize = function(...) { return (as.data.table(cm_backend_optimize(...))) }

# Load likelihood function for additional data
source(path("Analyses/2-suscsymp.R"));

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

expected_inf_dist = function(contact, pop, d_Ip, d_Is, d_Ia, f_Ip, f_Is, f_Ia, susc, symp)
{
    ngm = susc * t(t(contact) * (
        symp * (f_Ip * d_Ip + f_Is * d_Is) + 
        (1 - symp) * f_Ia * d_Ia)
    )
    dist = pop * abs(eigen(ngm)$vectors[,1])
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
likelihood = function(x, D, iter)
{
    ll_data = 0;
    fIa = D$f_Ia;

    if (D$f_Ia == "fit") {
        fIa = x[length(x)];
        x = x[1:(length(x) - 1)];
    }
    if (D$method != "multinomial") {
        sizeC = x[length(x)];
        x = x[1:(length(x) - 1)];
    }
    
    mult = x[length(x) - 1];
    sizeX = x[length(x)];
    x = x[1:(length(x) - 2)];
    
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
            symp = D$symp;
            path = x;
            path2 = NULL;
        } else if (D$vary == "symp") {
            susc = rep(1, n_age_groups);
            symp = c(rep(x, each = 2), rep(x[length(x)], n_age_groups - length(x)*2));
            path = x;
            path2 = NULL;
        } else if (D$vary == "both") {
            susc = c(rep(x[1:8], each = 2), rep(x[8], n_age_groups - 16));
            symp = c(rep(x[9:16], each = 2), rep(x[16], n_age_groups - 16));
            path = x[1:8];
            path2 = x[9:16];
        } else {
            stop("likelihood: vary must be 'susc' or 'symp' or 'both'.");
        }
        
        # expected case distribution given model parameters
        dist = expected_case_dist(d$mat, d$pop, 
            D$d_Ip, D$d_Is, D$d_Ia, 
            D$f_Ip, D$f_Is, fIa,
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
            ll_data = ll_data + sum(ddirmnom(counts, sum(counts), sizeC * dist, log = T)) * d$weight;
        }
    }
    
    # Fit to additional data
    ll_extra = likelihood_extra(x_susc = x[1:8], x_symp = x[9:16], x_subclin_mult = mult, x_size = sizeX, D$tag)
    
    return (ll_data + ll_extra)
}

# Do fitting
do_fit = function(age_dist, varying, locations, method, tag, fIa, burn_in = 3000, iter = 2000)
{
    # Get data for fitting
    fit = NULL;
    d = list(
        vary = varying,
        symp = 0.5,
        data = lapply(locations, make_data),
        d_Ip = 2.1,
        d_Is = 2.9,
        d_Ia = 5,
        f_Ip = 1,
        f_Is = 1,
        f_Ia = fIa,
        g_sd = 0.05,
        method = method,
        tag = tag
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
    if (varying == "both") {
        priors = c(make_priors(8, "B 1.5 1.5 I 0.5 0.505", "u_"), make_priors(8, "B 1.5 1.5 I 0.5 0.505", "y_"));
    } else {
        priors = make_priors(8, "B 2 2");
    }
    
    priors = c(priors, mult = "N 5 3 T 1 10");
    priors = c(priors, sizeX = "N 0 100 T 0 1000");
    
    if (method != "multinomial") {
        priors = c(priors, sizeC = "N 0 100 T 0 1000");
    }

    if (fIa == "fit") {
        priors = c(priors, fIa = "U 0 1");
    }
    
    fit = RMCMC(likelihood, d, priors, 0, burn_in, length(priors) * 2, iter, T, F, F, 1);

    return (fit)
}

all_locations = age_dist[, unique(location)];

argv = commandArgs();
argc = length(argv);
fIa = argv[argc];
if (fIa %like% "^[.0-9]+$") {
    fIa = as.numeric(fIa);
}
tag_ind = argv[argc-2];
tag_all = argv[argc-1];

if (tag_ind != "NULL") {

f_ind = NULL;
for (loc in "Hubei_CCDC") {#-#all_locations) {
    cat(paste0(loc, "...\n"));
    f = do_fit(age_dist, "both", loc, "multinomial", tag_ind, fIa);#-#
    f_ind = rbind(f_ind, cbind(f, location = loc))
}
qsave(f_ind, path(paste0("Analyses/2-HUBEI-linelist_mn_both_fit_ind_fIa", fIa, "-", tag_ind, ".qs")), "balanced");

# melted = melt(f_ind, id.vars = c("trial", "lp", "chain", "ll", "location"))
# theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))
# p1 = ggplot(melted[variable %like% "u_"]) + geom_violin(aes(x = variable, y = value)) + ylim(0, 1) + labs(y = NULL, x = NULL) + facet_wrap(~location, ncol = 1)
# p2 = ggplot(melted[variable %like% "y_"]) + geom_violin(aes(x = variable, y = value)) + ylim(0, 1) + labs(y = NULL, x = NULL) + facet_wrap(~location, ncol = 1)
# p3 = ggplot(melted[variable %like% "mult"]) + geom_violin(aes(x = variable, y = value)) + ylim(0, 10) + labs(y = NULL, x = NULL) + facet_wrap(~location, ncol = 1)
# #-#p4 = ggplot(melted[variable %like% "size"]) + geom_violin(aes(x = variable, y = value)) + ylim(0, 100) + labs(y = NULL, x = NULL) + facet_wrap(~location, ncol = 1)
# #-#p = plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(1, 1, 0.4, 0.4))
# p = plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, 1, 0.4))#-#
# ggsave(path(paste0("Analyses/2-mn-newfit-", tag_ind, ".pdf")), p, width = 15, height = 40, units = "cm", useDingbats = F)
}

if (tag_all != "NULL") {
# fit both
f_both = do_fit(age_dist, "both", all_locations, "dirichlet-multinomial", tag_all, fIa, 5000, 5000)
f_both = f_both[trial %% 10 == 0]
qsave(f_both, path(paste0("Analyses/2-linelist_both_fit_fIa", fIa, "-", tag_all, ".qs")), "balanced");
melted = melt(f_both, id.vars = c("trial", "lp", "chain", "ll"))

p1 = ggplot(melted[variable %like% "u_"]) + geom_violin(aes(x = variable, y = value), scale = "width") + ylim(0, 1) + labs(y = NULL, x = NULL);
p2 = ggplot(melted[variable %like% "y_"]) + geom_violin(aes(x = variable, y = value), scale = "width") + ylim(0, 1) + labs(y = NULL, x = NULL);
p3 = ggplot(melted[variable %like% "mult"]) + geom_violin(aes(x = variable, y = value), scale = "width") + ylim(0, 10) + labs(y = NULL, x = NULL);
p4 = ggplot(melted[variable %like% "size"]) + geom_violin(aes(x = variable, y = value), scale = "width") + ylim(0, 300) + labs(y = NULL, x = NULL);

p = plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(1, 1, 0.4, 0.4))
ggsave(path(paste0("Analyses/2-newfit-", tag_all, ".pdf")), p, width = 25, height = 10, units = "cm", useDingbats = F)
}

if (tag_ind != "NULL") {

f_ind = qread(paste0("~/Dropbox/nCoV/Analyses/2-linelist_mn_both_fit_ind_fIa0.5-", tag_ind, ".qs"));

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
    } else if (D$vary == "both") {
        susc = c(rep(x[1:8], each = 2), rep(x[8], n_age_groups - 16));
        symp = c(rep(x[9:16], each = 2), rep(x[16], n_age_groups - 16));
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
    if (length(dist) > length(x) / 2) {
        dist = c(dist[1:((length(x) / 2)-1)], sum(dist[(length(x) / 2):length(dist)]));
    }
    
    return (dist)
}

results = NULL;
for (loc in all_locations)
{
    d = list(
        vary = "both",
        symp = 0.5,
        data = lapply(loc, make_data),
        d_Ip = 2.1,#
        d_Is = 2.9,#
        d_Ia = 5,#
        f_Ip = 1,
        f_Is = 1,
        f_Ia = fIa,
        g_sd = 0.15
    );

    draws = NULL
    data = f_ind[location == loc];
    data = data[sample.int(.N, 100)];
    for (i in 1:nrow(data)) {
        cd = get_expected(unlist(unname(data[i, 5:20])), d);
        names(cd) = names(data[i, 5:12]);
        draws = rbind(draws, as.data.table(as.list(cd)));
    }
    
    results = rbind(results,
        cbind(draws, location = loc));
};

qsave(results, path(paste0("Analyses/2-linelist_both_expected_fIa", fIa, "-", tag_ind, ".qs")), "balanced")
}