// mcmc.cpp

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include "nloptrAPI.h"
#include "randomizer.h"
#include "distribution.h"
#include "mcmc.h"
using namespace std;

class Caller
{
public:
    Caller(Rcpp::Function ff, Rcpp::List ep)
     : f(ff), p(ep) {}

    double operator()(const vector<double>& theta, int& obs)
    {
        (void)obs;
        return Rcpp::as<double>(f(theta, p));
    }

private:
    Rcpp::Function f;
    Rcpp::List p;
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
Rcpp::DataFrame MCMC(Rcpp::Function likelihood, Rcpp::List extra_params,
    unsigned int burn_in, unsigned int iterations, unsigned int n_chains, 
    Rcpp::List params_priors)
{
    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    Reporter rep(iterations, n_chains, params);

    DEMCMC_Priors<int>(rand, lik, rep, burn_in, iterations, n_chains, priors, true, params);

    return rep.Contents();
}

// [[Rcpp::export]]
Rcpp::DataFrame MCMC_Init(Rcpp::Function likelihood, Rcpp::List extra_params,
    unsigned int burn_in, unsigned int iterations, unsigned int n_chains, 
    Rcpp::List params_priors, Rcpp::NumericMatrix initial)
{
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

    DEMCMC_Priors<int>(rand, lik, rep, burn_in, iterations, n_chains, priors, true, params,
        false, -1, init);

    return rep.Contents();
}

// [[Rcpp::export]]
Rcpp::List Optimize(Rcpp::Function likelihood, Rcpp::List extra_params,
    Rcpp::List params_priors, unsigned int global_iter = 1000, bool local = true)
{
    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    Reporter rep(1, 1, params);

    Optimize_Priors<int>(rand, lik, rep, priors, global_iter, local);

    return rep.ListContents();
}
