library(data.table)

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

gather(0.0, F)
gather(0.5, F)
gather(1.0, F)
gather(0.5, T)

