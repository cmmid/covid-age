// randomizer.h
// (C) 2013-2020 Nicholas G Davies

#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <random>
#include <vector>
#include <limits>
#include <gsl/gsl_rng.h>

class Randomizer
{
public:
    typedef std::mt19937_64 engine_type;
    engine_type engine;

    Randomizer();
    ~Randomizer();

    void Reset();

    double Uniform(double min = 0.0, double max = 1.0);
    double RoundedUniform(double min = 0.0, double max = 1.0, double shoulder = 0.01);
    double Normal(double mean = 0.0, double sd = 1.0);
    double Normal(double mean, double sd, double clamp);
    double Cauchy(double x0 = 0.0, double gamma = 1.0);
    double LogNormal(double mean = 0.0, double sd = 1.0);
    double Exponential(double rate = 1.0);
    double Gamma(double alpha, double beta);
    double Beta(double alpha, double beta);
    unsigned int Discrete(unsigned int size);
    int Discrete(int min, int max);
    int Discrete(std::vector<unsigned int>& cumulative_weights);
    int Discrete(std::vector<double>& cumulative_weights);
    void Multinomial(unsigned int N, std::vector<double>& p, std::vector<unsigned int>& n);
    int FlipCoin();
    void Pick(int min, int max, int n, std::vector<int>& picks);
    bool Bernoulli(double p);
    unsigned int Binomial(unsigned int n, double p);
    unsigned int BetaBinomial(unsigned int n, double p, double a_plus_b);
    int Poisson(double mean);
    int Geometric(double p);
    int NonzeroPoisson(double mean);
    int FoundressDual(double a, int n_max);
    // Approximate. Accurate within ~0.025% for a = 0, a = 1, and 0.01 < a < 0.99.
    int FoundressPoisson(double a, int n_max);
    int Round(double x);
    void SetEventRate(unsigned int handle, double p);
    bool Event(unsigned int handle);
    template <typename RandomAccessIterator>
    void Shuffle(RandomAccessIterator first, RandomAccessIterator last);

    void DiehardOutput(const char* filename);

    inline gsl_rng* GSL_RNG() { return r; }

private:
    double LambertW0(const double x);
    static double FoundressPoisson_LogLambda[256];

    std::seed_seq seed;
    gsl_rng * r; // for GSL random routines...

    engine_type::result_type fast_bits;
    int fast_shift;

    std::vector<std::geometric_distribution<unsigned int>> event_distributions;
    std::vector<unsigned int> steps_to_next_event;
};

template <typename RandomAccessIterator>
void Randomizer::Shuffle(RandomAccessIterator first, RandomAccessIterator last)
{
    auto n = last - first;
    for (size_t i = 0, size = last - first; i < n; ++i)
        std::swap(*(first + i), *(first + i + Discrete(size - i)));
}


#endif