`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/nick/Documents/ncov_age/covidm/build/sourceCpp-x86_64-apple-darwin15.6.0-1.0.4/sourcecpp_5ac168489ae8/sourceCpp_2.so')

cm_backend_simulate <- Rcpp:::sourceCppFunction(function(parameters, n_run = 1L, seed = 0L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cm_backend_simulate')
cm_evaluate_distribution <- Rcpp:::sourceCppFunction(function(dist_code, steps = 101L, xmin = 0, xmax = -1) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cm_evaluate_distribution')

rm(`.sourceCpp_1_DLLInfo`)
