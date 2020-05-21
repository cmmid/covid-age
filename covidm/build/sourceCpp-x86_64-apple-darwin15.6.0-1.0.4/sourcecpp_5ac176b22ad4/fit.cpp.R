`.sourceCpp_3_DLLInfo` <- dyn.load('/Users/nick/Documents/ncov_age/covidm/build/sourceCpp-x86_64-apple-darwin15.6.0-1.0.4/sourcecpp_5ac176b22ad4/sourceCpp_4.so')

cm_backend_mcmc <- Rcpp:::sourceCppFunction(function(likelihood, extra_params, params_priors, seed, burn_in, n_chains, iterations, verbose, reeval_likelihood, in_parallel, n_threads) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_cm_backend_mcmc')
cm_backend_mcmc_init <- Rcpp:::sourceCppFunction(function(likelihood, extra_params, params_priors, initial, seed, burn_in, n_chains, iterations, verbose, reeval_likelihood, in_parallel, n_threads) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_cm_backend_mcmc_init')
cm_backend_optimize <- Rcpp:::sourceCppFunction(function(likelihood, extra_params, params_priors, seed, global, global_algorithm, global_maxeval, global_ftol_abs, local, local_algorithm, local_maxeval, local_ftol_abs, verbose) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_cm_backend_optimize')
cm_backend_prior_sample <- Rcpp:::sourceCppFunction(function(params_priors) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_cm_backend_prior_sample')

rm(`.sourceCpp_3_DLLInfo`)
