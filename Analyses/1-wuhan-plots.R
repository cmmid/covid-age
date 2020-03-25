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

path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
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



calc_eigen = function(z, setter) {
    posterior_mean = colMeans(z[, 5:ncol(z)]);
    p = setter(posterior_mean);

    po = p$populations[[1]];
    d_Ip = sum(po$dist_Ip * seq(0, by = p$time_step, length.out = length(po$dist_Ip)));
    d_Is = sum(po$dist_Is * seq(0, by = p$time_step, length.out = length(po$dist_Is)));
    d_Ia = sum(po$dist_Ia * seq(0, by = p$time_step, length.out = length(po$dist_Ia)));
    
    ngm = po$susc * t(t(po$contact) * (
        po$symptom * (po$trans_Ip * d_Ip + po$trans_Is * d_Is) + 
        (1 - po$symptom) * po$trans_Ia * d_Ia)
    )
    
    eig = eigen(ngm)
    
    list(
        R0 = abs(eig$values[1]),
        stat = abs(eig$vectors[,1]),
        statc = abs(eig$vectors[,1]) * po$pop_size
    )
}

calc_R0_sub = function(p, contact, susc_frac = 1) {
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

calc_R0 = function(fit, n_samples) {
    s = sample.int(nrow(fit$posterior), n_samples);
    R0 = NULL;
    for (r in s) {
        # Get parameters
        x = unlist(fit$posterior[r,]);
        p = cm_translate_parameters(fit$parameters_func(fit$base_parameters, x));
        m = p$pop[[1]]$matrices;
        
        summat = function(m, con) { m[[1]] * con[[1]] + m[[2]] * con[[2]] + m[[3]] * con[[3]] + m[[4]] * con[[4]] }

        R0 = rbind(R0, data.table(t = 0, R0 = calc_R0_sub(p, summat(m, p$pop[[1]]$contact))));
        for (i in 1:length(p$pop[[1]]$schedule)) {
            R0 = rbind(R0, data.table(t = p$pop[[1]]$schedule[[i]]$t, R0 = calc_R0_sub(p, summat(m, p$pop[[1]]$schedule[[i]]$contact))));
        }
    }
    
    R0 = R0[, as.list(c(hdi(R0, 0.95), hdi(R0, 0.5), mean(R0))), by = t];
    names(R0) = c("t0", "lower95", "upper95", "lower50", "upper50", "mean");
    R0[, t1 := c(t0[2:length(t0)], p$time1)];
    return (R0)
}

show_R0 = function(R0, seed_date, time_start = "2019-11-01") {
    rr = rbind(R0[, .(            t0, scenario, lower95, upper95, lower50, upper50, mean)],
               R0[, .(t0 = t1 - 1e-6, scenario, lower95, upper95, lower50, upper50, mean)]);
    setkey(rr, t0);
    rr[, t := ymd(time_start) + t0];
    
    R0labs = R0[, .(scenario = scenario, t = ymd(time_start) + (t0 + t1)/2, y = upper95 + 0.2, value = paste0(round(mean, 1), " (", round(lower95, 1), "-", round(upper95, 1), ")"))]
    
    ggplot(rr) + 
        geom_histogram(data = seed_date, aes(x = d, y = 6 * stat(ndensity)), binwidth = 1, fill = "#ffcccc") +
        geom_ribbon(aes(x = t, ymin = lower95, ymax = upper95, fill = scenario), alpha = 0.5) +
        geom_ribbon(aes(x = t, ymin = lower50, ymax = upper50, fill = scenario), alpha = 1) +
        geom_vline(data = R0[t0 > 0], aes(xintercept = t0 + ymd(time_start)), linetype = "33", colour = "grey", size = 0.25) +
        geom_text(data = R0labs, aes(x = t, y = y, label = value, colour = scenario), size = 1, hjust = 0.5, vjust = 0) +
        labs(x = "Date", y = NULL, title = expression(R[0])) +
        facet_wrap(~paste0("H[", scenario, "]"), strip.position = "left", ncol = 1, labeller = label_parsed) +
        sc_colours() +
        theme(legend.position = "none", strip.placement = "outside", 
            strip.background = element_blank(), strip.text.y.left = element_text(angle = 0))
}


show_hist = function(title, ages, data, bounds)
{
    d = NULL;
    for (i in 1:length(ages)) {
        d = rbind(d, data.table(age = rep(ages[i], length(data[[i]])), x = data[[i]]));
    }
    
    d[, age := factor(age, levels = ages)];
    
    ggplot(d) + 
        geom_density(aes(x = x, y = stat(ndensity), group = age, fill = age), colour = NA) + 
        labs(x = title, y = "") +
        xlim(bounds) +
        theme(legend.position = "none")
}

add_coss = function(p, max_age, colour, x0, y0, x1, y1, x2, y2, lt = "solid")
{
    len = max(length(x0), length(y0), length(x1), length(y1), length(x2), length(y2));
    x0 = rep_len(x0, len);
    y0 = rep_len(y0, len);
    x1 = rep_len(x1, len);
    y1 = rep_len(y1, len);
    x2 = rep_len(x2, len);
    y2 = rep_len(y2, len);
    
    n_samples = 1000;
    ages = c(0, seq(2.5, max_age - 2.5, by = 5), max_age);
    rows = sample.int(length(x0), n_samples, replace = T);
    mat = matrix(0, nrow = n_samples, ncol = length(ages));

    for (i in 1:n_samples)
    {
        r = rows[i];
        mat[i, ] = ifelse(ages < x1[r], 
            cm_interpolate_cos(ages, x0[r], y0[r], x1[r], y1[r]), 
            cm_interpolate_cos(ages, x1[r], y1[r], x2[r], y2[r]));
    }
    
    d = data.table(i = 1:length(ages), age = ages);
    d = d[, c(as.list(hdi(mat[, i], 0.95)), as.list(hdi(mat[, i], 0.5)), mean = mean(mat[, i]), age = age), by = i];
    
    names(d) = c("i", "lo95", "hi95", "lo50", "hi50", "mean", "age");
    
    p + 
        geom_ribbon(data = d, aes(x = age, ymin = lo95, ymax = hi95), fill = colour, alpha = 0.4) +
        geom_ribbon(data = d, aes(x = age, ymin = lo50, ymax = hi50), fill = colour, alpha = 0.4) +
        geom_line(data = d, aes(x = age, y = mean), colour = colour, linetype = lt)
}

calc_DIC = function(z, setter, likelihood)
{
    Ep = unname(unlist(colMeans(z[, 5:ncol(z)])));
    D_mean = -2 * likelihood(Ep, list(setter()));
    mean_D = z[, mean(-2 * lp)];
    
    2 * mean_D - D_mean;
}

redo_age = function(x, col = "age")
{
    x[, (col) := as.character(get(col))];
    x[, (col) := str_replace_all(get(col), " ", "")];
    x[, (col) := factor(get(col), levels = unique(get(col)))];
    return (x)
}

sp_comparison4 = function(results, time_start = "2019-12-01", time_end = "2040-01-01")
{
    inc = results[compartment == "cases_reported", .(n = sum(value)), by = .(t, run, scenario)];
    incHDI = inc[, c(as.list(hdi(n)), "mean" = mean(n), "median" = median(n)), by = .(t, scenario)];
    incHDI[, date := ymd(time_start) + t];
    incHDI[, scenario := factor(scenario, levels = unique(scenario))];
    
    plEpiMain = ggplot(incHDI[date <= ymd(time_end)]) +
        geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = scenario, group = scenario), alpha = 0.25) +
        geom_line(aes(x = date, y = mean, colour = scenario, group = scenario)) +
        geom_point(data = cases_wuhan1, aes(x = date, y = pmax(0.1, incidence)), size = 0.25) +
        geom_errorbarh(data = cases_wuhan2, aes(xmin = date0, xmax = date1, y = incidence/8, height = 0)) +
        sc_colours() +
        labs(x = NULL, y = "Incident clinical cases") +
        theme(legend.position = "none")

    incHDI[, lower := pmax(0.1, lower)]
    incHDI[, mean := pmax(0.1, mean)]
    incHDI[, upper := pmax(0.1, upper)]
    plEpiInset = ggplot(incHDI[date <= ymd(time_end)]) +
        geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = scenario, group = scenario), alpha = 0.25) +
        geom_line(aes(x = date, y = mean, colour = scenario, group = scenario)) +
        geom_point(data = cases_wuhan1, aes(x = date, y = pmax(0.1, incidence)), size = 0.25) +
        geom_errorbarh(data = cases_wuhan2, aes(xmin = date0, xmax = date1, y = incidence/8, height = 0)) +
        sc_colours() +
        scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c("0", "1", "10", "100", "1000")) +
        labs(x = NULL, y = NULL) + theme(legend.position = "none")
    
    plEpi = ggdraw(plEpiMain) + draw_plot(plEpiInset, .1, .2, .6, .7)
    
    # show age dist from Li paper
    model_dist = results[, cm_case_distribution(.SD, time_start, "2019-12-08", c("2019-12-31", "2020-01-11", "2020-01-22"), c(0, 15, 45, 65, 100), "cases_reported"), by = .(run, scenario)];
    model_data = model_dist[, cm_mean_hdi(fcases), by = .(date, age, scenario)]
    model_data = rbind(model_data, cbind(age_dist_wuhan[date %in% unique(model_data$date)], scenario = "observed"))
    model_data[, scenario := factor(scenario, levels = c("observed", "1", "2", "3"))];
    
    model_data[date == "2019-12-08 - 2019-12-31", date := "Dec 8-31, 2019"];
    model_data[date == "2020-01-01 - 2020-01-11", date := "Jan 1-11, 2020"];
    model_data[date == "2020-01-12 - 2020-01-22", date := "Jan 12-22, 2020"];
    model_data = redo_age(model_data);

    plDist = ggplot(model_data) +
        geom_col(data = model_data[scenario != "observed"], 
            aes(x = age, y = 100 * upper, fill = scenario, group = scenario), colour = NA, alpha = 0.5,
            position = position_dodge(width = 0.8), width = 0.7) +
        geom_col(data = model_data[scenario != "observed"], 
            aes(x = age, y = 100 * lower, fill = scenario, group = scenario), colour = NA, 
            position = position_dodge(width = 0.8), width = 0.7) +
        geom_point(data = model_data[scenario != "observed"], 
            aes(x = age, y = 100 * mean, colour = scenario, group = scenario), 
            size = 0.25, position = position_dodge(width = 0.8), shape = 20) +
        geom_col(data = model_data[scenario == "observed"],
            aes(x = age, y = 100 * mean, fill = scenario, group = scenario, colour = scenario),
            size = 0.125) +
        geom_errorbar(data = model_data[scenario == "observed"],
            aes(x = age, ymin = 100 * lower, ymax = 100 * upper, group = scenario),
            size = 0.125, width = 0.2, colour = "black") +
        facet_wrap(~date, ncol = 3) +
        sc_colours("fill") + sc_outline_colours() +
        labs(x = "Age", y = "Clinical cases (%)") +
        theme(legend.position = "none", strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

    # again for China CDC
    model_dist = results[, cm_case_distribution(.SD, time_start, "2019-12-08", "2020-02-11", c(seq(0, 70, by = 10), 100), "cases_reported"), by = .(run, scenario)];
    model_data = model_dist[, cm_mean_hdi(fcases), by = .(date, age, scenario)]
    model_data = rbind(model_data, cbind(age_dist_wuhan[date %in% unique(model_data$date)], scenario = "observed"))
    model_data[, scenario := factor(scenario, levels = c("observed", "1", "2", "3"))];
    
    model_data[date == "2019-12-08 - 2020-02-11", date := "Dec 8, 2019 - Feb 11, 2020"];
    model_data = redo_age(model_data);

    plDist2 = ggplot(model_data) +
        geom_col(data = model_data[scenario != "observed"], 
            aes(x = age, y = 100 * upper, fill = scenario, group = scenario), colour = NA, alpha = 0.5,
            position = position_dodge(width = 0.8), width = 0.7) +
        geom_col(data = model_data[scenario != "observed"], 
            aes(x = age, y = 100 * lower, fill = scenario, group = scenario), colour = NA, 
            position = position_dodge(width = 0.8), width = 0.7) +
        geom_point(data = model_data[scenario != "observed"], 
            aes(x = age, y = 100 * mean, colour = scenario, group = scenario), 
            size = 0.25, position = position_dodge(width = 0.8), shape = 20) +
        geom_col(data = model_data[scenario == "observed"], 
            aes(x = age, y = 100 * mean, fill = scenario, group = scenario, colour = scenario), 
            size = 0.125, position = position_dodge(width = 0.9)) +
        geom_errorbar(data = model_data[scenario == "observed"], 
            aes(x = age, ymin = 100 * lower, ymax = 100 * upper, group = scenario), 
            size = 0.125, width = 0.2, position = position_dodge(width = 0.9), colour = "black") +
        facet_wrap(~date, ncol = 3) +
        sc_colours("fill") + sc_outline_colours() +
        labs(x = "Age", y = NULL) +
        theme(legend.position = "none", strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

    # again for subclinical cases
    model_dist = results[, cm_case_distribution(.SD, time_start, "2019-12-08", "2020-02-11", c(seq(0, 70, by = 10), 100), "subclinical"), by = .(run, scenario)];
    model_data = model_dist[, cm_mean_hdi(fcases), by = .(date, age, scenario)];
    model_data$scenario = as.character(model_data$scenario)
    model_data[, scenario := factor(scenario, levels = c("observed", "1", "2", "3"))];
    
    model_data[date == "2019-12-08 - 2020-02-11", date := "Dec 8, 2019 - Feb 11, 2020"];
    model_data = redo_age(model_data);

    plDistS = ggplot(model_data) +
        geom_col(aes(x = age, y = 100 * mean, fill = scenario, group = scenario, colour = scenario), size = 0.125, alpha = 0.5, position = position_dodge(width = 0.9)) +
        geom_errorbar(aes(x = age, ymin = 100 * lower, ymax = 100 * upper, group = scenario), size = 0.125, width = 0.2, position = position_dodge(width = 0.9), colour = "black") +
        facet_wrap(~date) +
        sc_colours("fill") + sc_outline_colours() +
        labs(x = "Age", y = "Subclinical cases (%)") +
        theme(legend.position = "none", strip.background = element_blank())

    plot_grid(plEpi, plot_grid(plDist, plDist2, nrow = 1, rel_widths = c(3, 2)), plDistS, labels = c("f", "g", "h"), label_size = 10, label_x = -0.025, nrow = 3, ncol = 1, rel_heights = c(2, 2, 2))
}

tiny_dist = function(d, title, col, tmax) {
    dd = data.table(y = d / max(d));
    dd[, t := (seq_along(d) - 1) * 0.25]
    mean = dd[, weighted.mean(t, y)]
    ggplot(dd) + 
        geom_ribbon(aes(x = t, ymin = 0, ymax = y), fill = col, colour = "black", size = 0.125) +
        annotate(geom = "text", x = tmax*2/3, y = 0.75, label = paste(round(mean, 1), "days"), size = 2) +
        labs(x = NULL, y = title) +
        coord_cartesian(xlim = c(0, tmax), expand = F) +
        theme(axis.title.y = element_text(angle = 0, hjust = 0, vjust = 0.5)) +
        scale_y_continuous(breaks = c(0, 1))
}

add_distr = function(d1, d2) {
    d = rep(0, length(d1));
    for (i1 in 1:length(d1)) {
        for (i2 in 1:length(d2)) {
            d[min(length(d), i1 + i2 - 1)] = d[min(length(d), i1 + i2 - 1)] + d1[i1] * d2[i2];
        }
    }
    return (d)
}

halve_distr = function(d1) {
    d = rep(0, length(d1));
    hl = (length(d1) - 1)/2
    d + c(d1[seq(1, length(d1), by = 2)], rep(0, hl)) + c(d1[seq(2, length(d1), by = 2)], rep(0, hl + 1));
}

show_posterior = function(fit)
{
    z = fit$posterior
    if (nrow(z) == 1) {
        z = rbind(z, z+1e-6, z-1e-6)
    }
    d = melt(z, id.vars = c("trial", "lp", "chain", "ll"));
    ggplot(d) + geom_density(aes(value, fill = variable)) + facet_wrap(~variable, scales = "free") + theme(legend.position = "none")
}

Col1 = "#ac2f2f";#"#b92d57";
Col2 = "#fb751e";#"#ffcc00";
Col3 = "#8ac3ff";#"#00b7ff";
sc_colours = function(a = c("colour", "fill")) scale_colour_manual(values = c("observed" = NA, "1" = Col1, "2" = Col2, "3" = Col3), aesthetics = a);
sc_outline_colours = function() scale_colour_manual(values = c("observed" = "#000000", "1" = Col1, "2" = Col2, "3" = Col3))
txt_theme = theme(plot.title = element_text(face = "plain", size = 7, hjust = 0.5))

###############
### FIGURES ###
###############

htitle = function(text, size = 6, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5) {
    ggdraw() + draw_label(text, x = x, y = y, size = size, hjust = hjust, vjust = vjust)
}


theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# FIG. 1 MODEL AND FITTING TO WUHAN

fig1 = function(fit_same, fit_susc, fit_symp, fileout, scf = 0.5)
{
    z_same = fit_same$posterior
    z_susc = fit_susc$posterior
    z_symp = fit_symp$posterior
    
    ppsusc = ggplot()
    ppsusc = add_coss(ppsusc, 90, Col1,
                    15, z_same$susc,
                    45, z_same$susc,
                    75, z_same$susc)
    ppsusc = add_coss(ppsusc, 90, Col2,
                    z_susc$age_y, z_susc$susc_y,
                    z_susc$age_m, z_susc$susc_m,
                    z_susc$age_o, z_susc$susc_o)
    ppsusc = add_coss(ppsusc, 90, Col3,
                    15, z_symp$susc,
                    45, z_symp$susc,
                    75, z_symp$susc)
    ppsusc = ppsusc + 
        ylim(0, 0.16) +
        labs(x = "Age", y = NULL, title = "Susceptibility") +
        txt_theme
    
    ppsymp = ggplot()
    ppsymp = add_coss(ppsymp, 90, Col1,
                    15, scf,
                    45, scf,
                    75, scf)
    ppsymp = add_coss(ppsymp, 90, Col2,
                    15, scf,
                    45, scf,
                    75, scf, "dashed")
    ppsymp = add_coss(ppsymp, 90, Col3,
                    z_symp$age_y, z_symp$symp_y,
                    z_symp$age_m, scf,
                    z_symp$age_o, z_symp$symp_o)
    ppsymp = ppsymp + 
        ylim(0, 1) +
        labs(x = "Age", y = NULL, title = "Clinical fraction") +
        txt_theme
    
    contact = rbind(
        data.table(scenario = "1", what = "q[H]", q = z_same$qH),
        data.table(scenario = "1", what = "q[L]", q = z_same$qL),
        data.table(scenario = "2", what = "q[H]", q = z_susc$qH),
        data.table(scenario = "2", what = "q[L]", q = z_susc$qL),
        data.table(scenario = "3", what = "q[H]", q = z_symp$qH),
        data.table(scenario = "3", what = "q[L]", q = z_symp$qL)
    )
    
    contact[, what := factor(what, levels = c("q[H]", "q[L]"))]
    
    ppcontact = ggplot(contact) +
        stat_ydensity(aes(what, q, fill = scenario), geom = "violin", colour = NA, kernel = "epanechnikov") +
        scale_x_discrete(NULL, limits = rev(levels(contact$what)), breaks = c("q[H]", "q[L]"), labels = parse(text = c("q[H]", "q[L]"))) +
        facet_wrap(~paste0("H[", scenario, "]"), strip.position = "left", ncol = 1, labeller = label_parsed) +
        sc_colours() +
        coord_flip() +
        labs(x = NULL, y = "Non-school contacts", title = "Contact multipliers") +
        theme(legend.position = "none", strip.placement = "outside", 
            strip.background = element_blank(), strip.text.y.left = element_text(angle = 0)) +
        txt_theme

    cat("Non-R0 done.\n")
    
    seed_date = rbind(
        data.table(scenario = "1", d = ymd(fit_same$base_parameters$date0) + fit_same$posterior$seed_start),
        data.table(scenario = "2", d = ymd(fit_susc$base_parameters$date0) + fit_susc$posterior$seed_start),
        data.table(scenario = "3", d = ymd(fit_symp$base_parameters$date0) + fit_symp$posterior$seed_start)
    );
    
    r0 = rbind(
        cbind(scenario = "1", calc_R0(fit_same, 100)),
        cbind(scenario = "2", calc_R0(fit_susc, 100)),
        cbind(scenario = "3", calc_R0(fit_symp, 100))
    )

    ppR0 = show_R0(r0, seed_date, time_start = fit_same$base_parameters$date0) + 
        ylim(0, 6) + txt_theme
    cat("R0 done.\n")

    pp = plot_grid(ppsusc, ppsymp,
        ppcontact, ppR0, nrow = 2, ncol = 2, align = "hv", axis = "bl", labels = c("b", "c", "d", "e"), label_size = 10)

    cat("Params done.\n")
    
    z = rbind(
        cbind(cm_sample_fit(fit_same, 50), scenario = 1),
        cbind(cm_sample_fit(fit_susc, 50), scenario = 2),
        cbind(cm_sample_fit(fit_symp, 50), scenario = 3)
    )
    
    pd = sp_comparison4(z, time_start = "2019-11-01", time_end = "2020-02-01")
    
    diagram = ggdraw() +
        draw_image(path("../model_diagram_small.png"));
    
    p1 = plot_grid(plot_grid(diagram, pp, nrow = 2, rel_heights = c(2, 4), labels = c("a", ""), label_size = 10), pd, nrow = 1, rel_widths = c(2, 3))

    cat("Epi curves done.\n")

    ggsave(path(paste0(fileout, ".pdf")), p1, width = 20, height = 10, units = "cm", useDingbats = F)
    ggsave(path(paste0(fileout, ".png")), p1, width = 20, height = 10, units = "cm")
}

tinydists = function(fit_same)
{
    pop1 = fit_same$base_parameters$pop[[1]];

    schematic = 
        plot_grid(
            plot_grid(
                diagram, 
                plot_grid(
                    tiny_dist(add_distr(pop1$dE, pop1$dIp), expression(P[I]), "#eeeeff", 21),
                    tiny_dist(add_distr(pop1$dE, halve_distr(0.5 * add_distr(pop1$dIp, pop1$dIs) + 0.5 * pop1$dIa)), expression(P[S]), "#ffeeff", 21), 
                ncol = 1),
            ncol = 2, rel_widths = c(3, 1)),
            plot_grid(
                tiny_dist(pop1$dE, expression(d[E]), "#ffeeee", 21),
                tiny_dist(pop1$dIp, expression(d[P]), "#ffffee", 21),
                tiny_dist(pop1$dIs, expression(d[C]), "#eeffee", 21),
                tiny_dist(pop1$dIa, expression(d[S]), "#eeffff", 21), ncol = 4),
            ncol = 1, rel_heights = c(2, 1)
        )
}


### PLOTTING

fit_same = cm_load(path("1-wuhan-fit-same-fIa-0.5.qs"))
fit_susc = cm_load(path("1-wuhan-fit-susc-fIa-0.5.qs"))
fit_symp = cm_load(path("1-wuhan-fit-symp-fIa-0.5.qs"))

fig1(fit_same, fit_susc, fit_symp, "../Submission/Fig1-test", 0.5);

cm_DIC(fit_same)
cm_DIC(fit_susc)
cm_DIC(fit_symp)

melt(fit_same$posterior, measure.vars = 5:ncol(fit_same$posterior))[, cm_mean_hdi(value), by = variable][
    , paste0(signif(mean,2), " (", signif(lower, 2), "-", signif(upper, 2), ")"), by = variable]

melt(fit_susc$posterior, measure.vars = 5:ncol(fit_susc$posterior))[, cm_mean_hdi(value), by = variable][
    , paste0(signif(mean,2), " (", signif(lower, 2), "-", signif(upper, 2), ")"), by = variable]

melt(fit_symp$posterior, measure.vars = 5:ncol(fit_symp$posterior))[, cm_mean_hdi(value), by = variable][
    , paste0(signif(mean,2), " (", signif(lower, 2), "-", signif(upper, 2), ")"), by = variable]


p = cm_plot_posterior(fit_same)
ggsave(path("../Submission/S1-wuhan-posterior-same.pdf"), p, width = 10, height = 10, units = "cm", useDingbats = F)
p = cm_plot_posterior(fit_susc)
ggsave(path("../Submission/S1-wuhan-posterior-susc.pdf"), p, width = 10, height = 10, units = "cm", useDingbats = F)
p = cm_plot_posterior(fit_symp)
ggsave(path("../Submission/S1-wuhan-posterior-symp.pdf"), p, width = 10, height = 10, units = "cm", useDingbats = F)




# POSTERIORS
cm_plot_posterior(cm_load(path("1-wuhan-fit-same-fIa-0.5.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-same-0.5.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
cm_plot_posterior(cm_load(path("1-wuhan-fit-susc-fIa-0.5.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-susc-0.5.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
cm_plot_posterior(cm_load(path("1-wuhan-fit-symp-fIa-0.5.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-symp-0.5.pdf", width = 12, height = 8, units = "cm", useDingbats = F)

cm_plot_posterior(cm_load(path("1-wuhan-fit-same-fIa-0.25.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-same-0.25.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
cm_plot_posterior(cm_load(path("1-wuhan-fit-susc-fIa-0.25.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-susc-0.25.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
cm_plot_posterior(cm_load(path("1-wuhan-fit-symp-fIa-0.25.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-symp-0.25.pdf", width = 12, height = 8, units = "cm", useDingbats = F)

cm_plot_posterior(cm_load(path("1-wuhan-fit-same-fIa-0.75.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-same-0.75.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
cm_plot_posterior(cm_load(path("1-wuhan-fit-susc-fIa-0.75.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-susc-0.75.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
cm_plot_posterior(cm_load(path("1-wuhan-fit-symp-fIa-0.75.qs")))
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/posterior-symp-0.75.pdf", width = 12, height = 8, units = "cm", useDingbats = F)
