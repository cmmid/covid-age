# 2-suscsymp.R: constraining susceptibility and clinical fraction curve with additional data.

# Load data

# Vo' matrix
m = cm_matrices[["Italy | Veneto"]]
matVo = m$home + m$work + m$other

# Vo' data
vo = data.table(read_excel("~/Dropbox/nCoV/vo.xlsx", 1));
vo[positive == T, sum(symptomatic == "yes")/.N, keyby = age_group]
vo[first_sampling == "Positive", .(symptomatic = sum(symptomatic == "yes"), .N), keyby = age_group]
vo[, symp1 := ifelse(symptomatic == "yes" & first_sampling == "Positive", T, F)]
vo[, fever_temperature := as.numeric(fever_temperature)]

# all italy data
italy_severity = fread(path("riccardo_symp_profile_it.csv"));
italy_severity = rbind(
    cbind(age = "0-9", round(italy_severity[1, 2:7] + italy_severity[2, 2:7] + italy_severity[3, 2:7]*(3/13))),
    cbind(age = "10-19", round(italy_severity[3, 2:7] * (10/13))),
    italy_severity[4:8],
    cbind(age = "70+", round(italy_severity[9, 2:7] + italy_severity[10, 2:7] + italy_severity[11, 2:7]))
)
italy_severity[, age := factor(age, levels = unique(age))]

# beta-binomial log density which allows for a range of possible successes
drbbinom = function(x_min, x_max, size, alpha, beta)
{
    # log(mean(dbbinom(x_min:x_max, size = size, alpha = alpha, beta = beta)))
    log(weighted.mean(dbbinom(x_min:x_max, size = size, alpha = alpha, beta = beta),
                      dbeta((0:(x_max - x_min) + 0.5)/(x_max - x_min + 1), 2, 2)))
}

# likelihood function
likelihood_extra = function(x_susc, x_symp, x_subclin_mult, x_size, tag)
{
    ll = 0
    
    # Fit to Gudbjartsson (Iceland) data on asymptomatic rate
    if (tag %like% "i") {
        # Coefficients here are the number of people in each age group who were tested.
        # It appears from Gudbjartsson et al. that not all test kits were returned, and this may be the number of people who 
        # were sent a test kit, as these data add up to 145 rather than 100. Assuming equal probability of returning tests 
        # regardless of age.
        iceland_analogue = (4 * x_symp[2] + 24 * x_symp[3] + 32 * x_symp[4] + 48 * x_symp[5] +
            27 * x_symp[6] + 9 * x_symp[7] + 1 * x_symp[8]) / 145;
        # Individuals with cough or fever: minimum 35, maximum 52
        ll = ll + drbbinom(35, 52, size = 100, alpha = iceland_analogue * x_size, beta = (1 - iceland_analogue) * x_size);
    }
    
    # Fit to Lavezzo (Vo) data on subclinical rate
    if (tag %like% "v") {
        # Coefficients here are the number of people in each age group who were tested; 80 in total.
        vo_analogue = (4 * x_symp[2] + 4 * x_symp[3] + 7 * x_symp[4] + 6 * x_symp[5] + 18 * x_symp[6] + 17 * x_symp[7] + 24 * x_symp[8]) / 80;
        # Individuals with cough or fever: minimum 35, maximum 45
        ll = ll + drbbinom(35, 45, size = 80, alpha = vo_analogue * x_size, beta = (1 - vo_analogue) * x_size);
    }
    
    # Fit to Nguyen (Ho Chi Minh City) data on subclinical rate
    if (tag %like% "h") {
        # Coefficients here are the number of infected people in each age group; we know there were 15 aged 16-29, and 15 aged 30-60
        hcmc_analogue = (15 * weighted.mean(x_symp[2:3], c(4, 10)) + 15 * mean(x_symp[4:6])) / 30;
        # Individuals with cough or fever: minimum 10, maximum 17
        ll = ll + drbbinom(10, 17, size = 30, alpha = hcmc_analogue * x_size, beta = (1 - hcmc_analogue) * x_size);
    }
    
    # Fit to Riccardo et al (Italy) data
    if (tag %like% "r") {
        italy_clinical = italy_severity[, (mild + severe + critical)];
        italy_subclinical = italy_severity[, (asymptomatic + paucisymptomatic)];
        italy_total = italy_severity[, total];
        italy_est_clinfrac = x_symp / ((1 - x_symp) / x_subclin_mult + x_symp);
        ll = ll + sum(dbbinom(italy_clinical, size = italy_total, alpha = x_size * italy_est_clinfrac, beta = x_size * (1 - italy_est_clinfrac), log = T));
    }
    
    # Fit to Bi et al data on susceptibility
    if (tag %like% "b") {
        bi_n = c(148, 85, 114, 268, 143, 110, 130, 72);
        bi_infected = c(11, 6, 7, 16, 7, 10, 20, 7);
        bi_overall_infected = sum(bi_infected) / sum(bi_n);
        ll = ll + sum(dbbinom(bi_infected, size = bi_n, alpha = bi_overall_infected * x_susc * x_size, beta = (1 - bi_overall_infected * x_susc) * x_size, log = T));
    }

    # Fit to Zhang et al data on susceptibility
    if (tag %like% "z") {
        zhang_n = c(305, 2210, 263);
        zhang_infected = c(8, 108, 17);
        zhang_overall_infected = sum(zhang_infected) / sum(zhang_n);
        # Hunan province age distribution by 10-year age group (millions):
        # 8.1 (0-9), 7.6, 10.9, 10.0, 11.7, 7.7, 5.2, 3.2, 1.1 (80+)
        x_susc3 = c(weighted.mean(x_susc[1:2], c(8.1, 3.8)), 
            weighted.mean(x_susc[2:7], c(3.8, 10.9, 10.0, 11.7, 7.7, 2.6)), 
            weighted.mean(x_susc[7:8], c(2.6, 4.4)));
        ll = ll + sum(dbbinom(zhang_infected, size = zhang_n, alpha = zhang_overall_infected * x_susc3 * x_size, beta = (1 - zhang_overall_infected * x_susc3) * x_size, log = T));
    }
    
    ll
}