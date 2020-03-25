# Load packages
library(data.table)
library(ggplot2)
library(stringr)
library(readxl)

# SHO LINE LIST: Provinces of China
sheets = excel_sheets("~/Dropbox/nCov-2019/data_sources/case_data/Linelist_20200211_3080.xlsx");
lls = NULL;
for (s in sheets) {
    l = read_excel("~/Dropbox/nCov-2019/data_sources/case_data/Linelist_20200211_3080.xlsx", s);
    lls = rbind(lls, l);
}
lls = data.table(lls);

do_lls = function(lls, prv) {
    counts = data.table(age_lower = seq(0, 90, by = 10));
    counts[, age_upper := age_lower + 10];
    counts[age_lower == 90, age_upper := 120];
    ll = lls[prv_en == prv & !is.na(age), .(age_lower = (age %/% 10) * 10)];
    ll[age_lower > 90, age_lower := 90];
    ll = cbind(location = prv, merge(counts, ll[, .(n = .N), by = age_lower], by = "age_lower", all = T))
}
 
# All provinces with >60 cases with defined ages 
lls_locs = c(
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

LL = NULL;
for (loc in lls_locs) {
    LL = rbind(LL, do_lls(lls, loc));
}

ggplot(LL[location == "Hubei"]) + geom_col(aes(age_lower, n))

# ALEX KOH: Singapore
sp = fread("~/Dropbox/nCoV/singapore_alexkoh_march_4.txt")
counts = data.table(age_lower = seq(0, 90, by = 10));
counts[, age_upper := age_lower + 10];
counts[age_lower == 90, age_upper := 120];
sp = sp[, .(age_lower = (age %/% 10) * 10)];
sp[age_lower > 90, age_lower := 90];
LL = rbind(LL, cbind(location = "Singapore", merge(counts, sp[, .(n = .N), by = age_lower], by = "age_lower", all = T)))

# SOUTH KOREA DATASET
sk = fread("~/Dropbox/nCoV/SKdataset/patient.csv")[!is.na(birth_year)];
sk[, age := 2020 - birth_year];
counts = data.table(age_lower = seq(0, 90, by = 10));
counts[, age_upper := age_lower + 10];
counts[age_lower == 90, age_upper := 120];
sk = sk[, .(age_lower = (age %/% 10) * 10)];
sk[age_lower > 90, age_lower := 90];
LL = rbind(LL, cbind(location = "South Korea", merge(counts, sk[, .(n = .N), by = age_lower], by = "age_lower", all = T)))

# OUTSIDE HUBEI LINELIST: South Korea, Japan
llo = fread("~/Dropbox/nCoV/ncov_outside_hubei.csv");
llo[, age0 := str_split_fixed(age, "-", 2)[,1]]
llo[, age1 := str_split_fixed(age, "-", 2)[,2]]
llo[age0 %like% "[0-9.]+", age2 := as.numeric(age0)]
llo[age0 %like% "[0-9.]+" & age1 %like% "[0-9.]+", age2 := runif(.N, as.numeric(age0), as.numeric(age1))]

# NOTE: Japan has 7 entries in this linelist with age recorded as 40-89.
# Based on inspection of the age distribution with these 7 cases removed,
# the most likely attribution of these cases is 2x 40, 2x 50, 2x 60, 1x 70.
# We do this rather than randomly assigning a value to avoid too much variability.
do_llo = function(llo, loc) {
    if (loc == "Japan") {
        llo = llo[age != "40-89"]
    }
    counts = data.table(age_lower = seq(0, 90, by = 10));
    counts[, age_upper := age_lower + 10];
    counts[age_lower == 90, age_upper := 120];
    ll = llo[country == loc & !is.na(age2), .(age_lower = (age2 %/% 10) * 10)];
    ll[age_lower > 90, age_lower := 90];
    if (loc == "Japan") {
        ll = rbind(ll, data.table(age_lower = c(40, 40, 60, 60, 70, 70, 80)))
    }
    ll = cbind(location = loc, merge(counts, ll[, .(n = .N), by = age_lower], by = "age_lower", all = T))
}

#LL = rbind(LL, do_llo(llo, "South Korea"))
LL = rbind(LL, do_llo(llo, "Japan"))


# CHINA CDC: WUHAN, HUBEI, REST OF CHINA
ccdc = data.table(read_excel("~/Dropbox/nCoV/ChinaCDC_agedist.xlsx", "Data"));
LL = rbind(LL, ccdc[location %in% c("Wuhan_CCDC", "Hubei_CCDC", "ChinaNonHubei_CCDC")]); ### TODO note that I feel this should be using Hubei outside Wuhan. Good to get a better estimate of the number of cases. Not sure where 20000 comes from.

# ITALY
LL = rbind(LL, fread("~/Dropbox/nCoV/italy-agedist-2020-03-20.csv"));

# ONTARIO
on = fread("~/Dropbox/nCoV/canada_processed.txt");
on = on[province == "Ontario" & age != "Not Reported"];
on[, age := as.numeric(str_replace(age, "^[<]?([0-9]+).*$", "\\1"))]; # sufficient for 10-year age bands. Assuming the one <18 is 10-19
on = on[, .(age_lower = (age %/% 10) * 10)];
counts = data.table(age_lower = seq(0, 90, by = 10));
counts[, age_upper := age_lower + 10];
counts[age_lower == 90, age_upper := 120];
LL = rbind(LL, cbind(location = "Ontario", merge(counts, on[, .(n = .N), by = age_lower], by = "age_lower", all = T)));

# replace NAs with 0s in counts
LL[is.na(n), n := 0]
fwrite(LL, "~/Dropbox/nCoV/Analyses/age.csv")

# END HERE

# # hubei linelist...? is biased away from children:
# llh = fread("~/Dropbox/nCoV/ncov_hubei.csv");
# llh[, age0 := str_split_fixed(age, "-", 2)[,1]]
# llh[, age1 := str_split_fixed(age, "-", 2)[,2]]
# llh[age0 %like% "[0-9.]+", age2 := as.numeric(age0)]
# llh[age0 %like% "[0-9.]+" & age1 %like% "[0-9.]+",
#     age2 := runif(.N, as.numeric(age0), as.numeric(age1))]
# ggplot(llh) + geom_histogram(aes(age2), binwidth = 10)
# 
# q = do_llo(llh, "China")
