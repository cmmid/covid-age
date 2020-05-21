// mcmc.cpp

// [[Rcpp::plugins(cpp11)]]
/* Rcpp::plugins(openmp) */
// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include "nloptrAPI.h"
#include "../randomizer/randomizer.h"
#include "../randomizer/distribution.h"
#include "mcmc.h"
using namespace std;

class Caller
{
public:
    Caller(Rcpp::Function ff, Rcpp::List ep)
     : f(ff), p(ep), i(0) {}

    double operator()(const vector<double>& theta, int& obs)
    {
        (void)obs;
        return Rcpp::as<double>(f(theta, p, ++i));
    }

private:
    Rcpp::Function f;
    Rcpp::List p;
    unsigned int i;
};

class Reporter
{
public:
    Reporter(unsigned int iterations, unsigned int chains, vector<string>& param_names)
     : n_chains(chains), n_samp(iterations * n_chains), 
       trial(n_samp, 0), chain(n_samp, 0), lp(n_samp, 0.0), ll(n_samp, 0.0), 
       theta(n_samp, param_names.size()), pnames(param_names)
    {
    }

    void operator()(int TR, double LP, double CH, double LL, vector<double> TH, int& OBS)
    {
        (void)OBS;

        unsigned int i = TR * n_chains + CH;

        if (TR < 0)
            i = 0;

        trial[i] = TR;
        chain[i] = CH;
        lp[i] = LP;
        ll[i] = LL;
        for (unsigned int j = 0; j < TH.size(); ++j)
            theta(i, j) = TH[j];
    }

    Rcpp::DataFrame Contents()
    {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("trial") = trial,
            Rcpp::Named("lp") = lp,
            Rcpp::Named("chain") = chain,
            Rcpp::Named("ll") = ll
        );

        for (unsigned int j = 0; j < pnames.size(); ++j)
            df.push_back(theta.column(j), pnames[j]);

        return Rcpp::DataFrame(df);
    }

    Rcpp::List ListContents()
    {
        Rcpp::List li = Rcpp::List::create(
            Rcpp::Named("trial") = trial[0],
            Rcpp::Named("lp") = lp[0],
            Rcpp::Named("chain") = chain[0],
            Rcpp::Named("ll") = ll[0]
        );

        for (unsigned int j = 0; j < pnames.size(); ++j)
            li.push_back(theta(j, 0), pnames[j]);

        return Rcpp::List(li);
    }

private:
    unsigned int n_chains, n_samp;
    Rcpp::IntegerVector trial, chain;
    Rcpp::NumericVector lp, ll;
    Rcpp::NumericMatrix theta;
    vector<string> pnames;
};

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_mcmc(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors,
    int seed, unsigned int burn_in, unsigned int n_chains, unsigned int iterations, bool verbose, 
    bool reeval_likelihood, bool in_parallel, int n_threads)
{
    /// TODO implement seed
    (void) seed;

    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    Reporter rep(iterations, n_chains, params);

    DEMCMC_Priors<int>(rand, lik, rep, burn_in, iterations, n_chains, priors, verbose, params, 
        reeval_likelihood, in_parallel, n_threads, false);

    return rep.Contents();
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_mcmc_init(Rcpp::Function likelihood, Rcpp::List extra_params,
    Rcpp::List params_priors, Rcpp::NumericMatrix initial,
    int seed, unsigned int burn_in, unsigned int n_chains, unsigned int iterations, bool verbose, 
    bool reeval_likelihood, bool in_parallel, int n_threads)
{
    /// TODO implement seed
    (void) seed;

    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    vector<vector<double>> init;

    if (initial.nrow() > 0)
    {
        for (unsigned int c = 0; c < n_chains; ++c)
        {
            if (initial.nrow() <= c)
                throw logic_error("initial list must be large enough to supply all chains");
            Rcpp::NumericVector nv = initial.row(c); 
            vector<double> v = Rcpp::as<vector<double>>(nv);
            if (v.size() != params_priors.size())
                throw logic_error("elements of initial list must be same size as number of parameters");
            init.push_back(v);
        }
    }

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    Reporter rep(iterations, n_chains, params);

    DEMCMC_Priors<int>(rand, lik, rep, burn_in, iterations, n_chains, priors, verbose, params, 
        reeval_likelihood, in_parallel, n_threads, true, -1, init);

    return rep.Contents();
}

// [[Rcpp::export]]
Rcpp::List cm_backend_optimize(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors,
    int seed, bool global, string global_algorithm, unsigned int global_maxeval, double global_ftol_abs,
    bool local, string local_algorithm, unsigned int local_maxeval, double local_ftol_abs, bool verbose)
{
    (void) seed;

    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    Reporter rep(1, 1, params);

    Optimize_Priors<int>(rand, lik, rep, priors,
        global, global_algorithm, global_maxeval, global_ftol_abs,
        local, local_algorithm, local_maxeval, local_ftol_abs, verbose);

    return rep.ListContents();
}

// [[Rcpp::export]]
Rcpp::NumericVector cm_backend_prior_sample(Rcpp::List params_priors)
{
    Randomizer rand;
    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    Rcpp::NumericVector prior_sample;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
    {
        Distribution dist(Rcpp::as<string>(params_priors[i]));
        prior_sample.push_back(dist.Random(rand), params[i]);
    }

    return prior_sample;
}


#include <Rcpp.h>
// cm_backend_mcmc
Rcpp::DataFrame cm_backend_mcmc(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors, int seed, unsigned int burn_in, unsigned int n_chains, unsigned int iterations, bool verbose, bool reeval_likelihood, bool in_parallel, int n_threads);
RcppExport SEXP sourceCpp_3_cm_backend_mcmc(SEXP likelihoodSEXP, SEXP extra_paramsSEXP, SEXP params_priorsSEXP, SEXP seedSEXP, SEXP burn_inSEXP, SEXP n_chainsSEXP, SEXP iterationsSEXP, SEXP verboseSEXP, SEXP reeval_likelihoodSEXP, SEXP in_parallelSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type likelihood(likelihoodSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type extra_params(extra_paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params_priors(params_priorsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_chains(n_chainsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type reeval_likelihood(reeval_likelihoodSEXP);
    Rcpp::traits::input_parameter< bool >::type in_parallel(in_parallelSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cm_backend_mcmc(likelihood, extra_params, params_priors, seed, burn_in, n_chains, iterations, verbose, reeval_likelihood, in_parallel, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// cm_backend_mcmc_init
Rcpp::DataFrame cm_backend_mcmc_init(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors, Rcpp::NumericMatrix initial, int seed, unsigned int burn_in, unsigned int n_chains, unsigned int iterations, bool verbose, bool reeval_likelihood, bool in_parallel, int n_threads);
RcppExport SEXP sourceCpp_3_cm_backend_mcmc_init(SEXP likelihoodSEXP, SEXP extra_paramsSEXP, SEXP params_priorsSEXP, SEXP initialSEXP, SEXP seedSEXP, SEXP burn_inSEXP, SEXP n_chainsSEXP, SEXP iterationsSEXP, SEXP verboseSEXP, SEXP reeval_likelihoodSEXP, SEXP in_parallelSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type likelihood(likelihoodSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type extra_params(extra_paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params_priors(params_priorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_chains(n_chainsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type reeval_likelihood(reeval_likelihoodSEXP);
    Rcpp::traits::input_parameter< bool >::type in_parallel(in_parallelSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cm_backend_mcmc_init(likelihood, extra_params, params_priors, initial, seed, burn_in, n_chains, iterations, verbose, reeval_likelihood, in_parallel, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// cm_backend_optimize
Rcpp::List cm_backend_optimize(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors, int seed, bool global, string global_algorithm, unsigned int global_maxeval, double global_ftol_abs, bool local, string local_algorithm, unsigned int local_maxeval, double local_ftol_abs, bool verbose);
RcppExport SEXP sourceCpp_3_cm_backend_optimize(SEXP likelihoodSEXP, SEXP extra_paramsSEXP, SEXP params_priorsSEXP, SEXP seedSEXP, SEXP globalSEXP, SEXP global_algorithmSEXP, SEXP global_maxevalSEXP, SEXP global_ftol_absSEXP, SEXP localSEXP, SEXP local_algorithmSEXP, SEXP local_maxevalSEXP, SEXP local_ftol_absSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type likelihood(likelihoodSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type extra_params(extra_paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params_priors(params_priorsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type global(globalSEXP);
    Rcpp::traits::input_parameter< string >::type global_algorithm(global_algorithmSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type global_maxeval(global_maxevalSEXP);
    Rcpp::traits::input_parameter< double >::type global_ftol_abs(global_ftol_absSEXP);
    Rcpp::traits::input_parameter< bool >::type local(localSEXP);
    Rcpp::traits::input_parameter< string >::type local_algorithm(local_algorithmSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type local_maxeval(local_maxevalSEXP);
    Rcpp::traits::input_parameter< double >::type local_ftol_abs(local_ftol_absSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cm_backend_optimize(likelihood, extra_params, params_priors, seed, global, global_algorithm, global_maxeval, global_ftol_abs, local, local_algorithm, local_maxeval, local_ftol_abs, verbose));
    return rcpp_result_gen;
END_RCPP
}
// cm_backend_prior_sample
Rcpp::NumericVector cm_backend_prior_sample(Rcpp::List params_priors);
RcppExport SEXP sourceCpp_3_cm_backend_prior_sample(SEXP params_priorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type params_priors(params_priorsSEXP);
    rcpp_result_gen = Rcpp::wrap(cm_backend_prior_sample(params_priors));
    return rcpp_result_gen;
END_RCPP
}
