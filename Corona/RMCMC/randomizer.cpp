// randomizer.cpp
// (C) 2013-2018 Nicholas G Davies
// Lambert W function implementation by Darko Veberic based on work by Toshio Fukushima (see below)

#include "randomizer.h"
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>


Randomizer::Randomizer()
 : seed({ 120368978, 37761590, 364135833, 399444367, 298076336 }), r(gsl_rng_alloc(gsl_rng_mt19937))
{
    Reset();
}

Randomizer::~Randomizer()
{
    gsl_rng_free(r);
}

void Randomizer::Reset()
{
    gsl_rng_set(r, 0); // default seed
    engine.seed(seed);
    fast_bits = engine();
    fast_shift = 0;
}

double Randomizer::Uniform(double min, double max)
{
    if (min == max)
        return min;
    std::uniform_real_distribution<double> dist(min, max);
    return dist(engine);
}

double Randomizer::RoundedUniform(double min, double max, double shoulder)
{
    if (min >= max)
        return min;
    double z = Uniform();
    double sd = shoulder * (max - min) / ((1 - shoulder) * 2.50662827463);
    if (z < shoulder / 2)
        return min - abs(Normal(0, sd));
    else if (z < shoulder)
        return max + abs(Normal(0, sd));
    else
        return Uniform(min, max);
}

double Randomizer::Normal(double mean, double sd)
{
    std::normal_distribution<double> dist(mean, sd);
    return dist(engine);
}

double Randomizer::Normal(double mean, double sd, double clamp)
{
    double n;
    do n = Normal(mean, sd); while (std::fabs(n - mean) > clamp);
    return n;
}

double Randomizer::LogNormal(double mean, double sd)
{
    std::lognormal_distribution<double> dist(mean, sd);
    return dist(engine);
}

double Randomizer::Cauchy(double x0, double gamma)
{
    std::cauchy_distribution<double> dist(x0, gamma);
    return dist(engine);
}

double Randomizer::Exponential(double rate)
{
    std::exponential_distribution<double> dist(rate);
    return dist(engine);
}

double Randomizer::Gamma(double alpha, double beta)
{
    std::gamma_distribution<double> dist(alpha, beta);
    return dist(engine);
}

double Randomizer::Beta(double alpha, double beta)
{
    double x = Gamma(alpha, 1.);
    double y = Gamma(beta, 1.);
    return x / (x + y);
}

unsigned int Randomizer::Discrete(unsigned int size)
{
    std::uniform_int_distribution<unsigned int> dist(0, size - 1);
    return dist(engine);
}

int Randomizer::Discrete(int min, int max)
{
    std::uniform_int_distribution<int> dist(min, max);
    return dist(engine);
}

int Randomizer::Discrete(std::vector<unsigned int>& cumulative_weights)
{
    if (cumulative_weights.back() == 0)
        return Discrete(cumulative_weights.size());
    return std::distance(cumulative_weights.begin(),
                std::lower_bound(cumulative_weights.begin(),
                    cumulative_weights.end(),
                    Discrete(cumulative_weights.back())));
}

int Randomizer::Discrete(std::vector<double>& cumulative_weights)
{
    if (cumulative_weights.back() == 0)
        return Discrete(cumulative_weights.size());
    return std::distance(cumulative_weights.begin(),
                std::lower_bound(cumulative_weights.begin(),
                    cumulative_weights.end(),
                    Uniform(0.0, cumulative_weights.back())));
}

void Randomizer::Multinomial(unsigned int N, std::vector<double>& p, std::vector<unsigned int>& n)
{
    gsl_ran_multinomial(r, p.size(), N, &p[0], &n[0]);
}

// Returns either 0 or 1 with equal probability.
int Randomizer::FlipCoin()
{
    if (fast_shift == std::numeric_limits<engine_type::result_type>::digits)
    {
        fast_bits = engine();
        fast_shift = 0;
    }
    return (fast_bits >> fast_shift++) & 1;
}


void Randomizer::Pick(int min, int max, int n, std::vector<int>& picks)
{
    picks.clear();

    for (int j = max - n + 1; j <= max; ++j)
    {
        int t = Discrete(min, j);
        auto ins = lower_bound(picks.begin(), picks.end(), t);
        if (ins == picks.end() || *ins != t)
            picks.insert(ins, t);
        else
            picks.push_back(j);
    }
}

bool Randomizer::Bernoulli(double p)
{
    if (p <= 0) return false;
    if (p >= 1) return true;
    std::bernoulli_distribution dist(p);
    return dist(engine);
}

unsigned int Randomizer::Binomial(unsigned int n, double p)
{
    if (p <= 0) return 0;
    return gsl_ran_binomial(r, p, n);
    //std::binomial_distribution<int> dist(n, p);
    //return dist(engine);
}

unsigned int Randomizer::BetaBinomial(unsigned int n, double p, double a_plus_b)
{
    if (a_plus_b > 0)
        p = gsl_ran_beta(r, a_plus_b * p, a_plus_b * (1 - p));
    return gsl_ran_binomial(r, p, n);
}

int Randomizer::Poisson(double mean)
{
    if (mean <= 0) return 0;
    std::poisson_distribution<int> dist(mean);
    return dist(engine);
}

int Randomizer::Geometric(double p)
{
    if (p <= 0) return 0;
    std::geometric_distribution<int> dist(p);
    return dist(engine);
}

int Randomizer::NonzeroPoisson(double mean)
{
    static double last_mean = -1, last_lambda = -1, last_t = -1;

    if (mean <= 1) return 1;

    // Solve for t and lambda to be used for zero-truncated Poisson distribution (cached)
    double t, lambda;
    if (mean == last_mean)
    {
        lambda = last_lambda;
        t = last_t;
    }
    else
    {
        lambda = mean + LambertW0(-mean * std::exp(-mean));
        t = std::exp(-lambda) / (1. - std::exp(-lambda)) * lambda;
        last_mean = mean;
        last_lambda = lambda;
        last_t = t;
    }

    // Generate a random variable from a zero-truncated Poisson distribution.
    // Follows Borje, Gio. "Zero-Truncated Poisson Distribution Sampling Algorithm".
    // http://giocc.com/zero_truncated_poisson_sampling_algorithm.html
    int k = 1;
    double s = t, u = Uniform();

    while (s < u)
    {
        ++k;
        t *= lambda / k;
        s += t;
    }

    return k;
}

int Randomizer::FoundressDual(double a, int n_max)
{
    if (a <= 0)
        return n_max;

    int n = std::floor(1 / a);
    double p = n * (a * (1 + n) - 1);
    int r = Bernoulli(p) ? n : n + 1;
    return std::min(n_max, r);
}

int Randomizer::FoundressPoisson(double a, int n_max)
{
    static double last_a = -1, last_lambda = -1, last_t = -1;

    if (a <= 0)
        return n_max;
    if (a >= 1)
        return 1;

    // Solve for t and lambda to be used for zero-truncated Poisson distribution, cached
    double t, lambda;
    if (a == last_a)
    {
        lambda = last_lambda;
        t = last_t;
    }
    else
    {
        double ai = 256 * a;

        auto quadratic_interpolate = [](double alpha, double beta, double gamma, double x)
        {
            double c = alpha;
            double a = (gamma - 2 * beta + alpha) / 2;
            double b = beta - alpha - a;
            return a * x * x + b * x + c;
        };

        if (ai < 16) // approximation based on eyeballing . . .
            lambda = 1 / a + pow(1.00393719292054584, 256 * a);
        else // approximation based on quadratic interpolation of log-lambda table . . .
        {
            int i = std::min((int)ai - 1, 253);
            lambda = exp(quadratic_interpolate(FoundressPoisson_LogLambda[i],
                                               FoundressPoisson_LogLambda[i + 1],
                                               FoundressPoisson_LogLambda[i + 2], ai - i - 1));
        }

        t = std::exp(-lambda) / (1. - std::exp(-lambda)) * lambda;

        last_a = a;
        last_lambda = lambda;
        last_t = t;
    }

    // Generate a random variable from a zero-truncated Poisson distribution.
    // Either does rejection sampling, or follows Gio Borje, "Zero-Truncated Poisson Distribution Sampling Algorithm".
    // http://giocc.com/zero_truncated_poisson_sampling_algorithm.html
    if (lambda > 16)
    {
        int ans;
        while (!(ans = Poisson(lambda)))
            ;
        return ans;
    }

    int k = 1;
    double s = t, u = Uniform();

    while (s < u)
    {
        ++k;
        t *= lambda / k;
        s += t;
    }

    return k;
}

int Randomizer::Round(double x)
{
    int sign = x < 0 ? -1 : 1;
    double intpart, fracpart;
    fracpart = std::modf(std::fabs(x), &intpart);
    return sign * (intpart + Bernoulli(fracpart));
}

// For a given positive handle, sets a probability p of trial success.
// Then the method Event(handle) below will return true with probability p.
void Randomizer::SetEventRate(unsigned int handle, double p)
{
    if (handle >= event_distributions.size())
    {
        event_distributions.resize(handle + 1);
        steps_to_next_event.resize(handle + 1);
    }

    event_distributions[handle] = std::geometric_distribution<unsigned int>(p);
    steps_to_next_event[handle] = event_distributions[handle](engine);
}

bool Randomizer::Event(unsigned int handle)
{
    if (event_distributions[handle].p() > 0 && steps_to_next_event[handle]-- == 0)
    {
        steps_to_next_event[handle] = event_distributions[handle](engine);
        return true;
    }
    return false;
}


void Randomizer::DiehardOutput(const char* filename)
{
    std::ofstream out(filename);
    for (int i = 0; i < 12000000;)
    {
        auto k = engine();
        out.write(reinterpret_cast<char*>(&k), sizeof(k));
        i += sizeof(k);
    }
}

// Adapted (lightly) from the implementation of the Lambert W function, 0 branch, by Darko Veberic,
// based on the method of Toshio Fukushima; see https://github.com/DarkoVeberic/LambertW
double LambertWSeries(const double p)
{
  static const double q[] = {
    -1,
    +1,
    -0.333333333333333333,
    +0.152777777777777778,
    -0.0796296296296296296,
    +0.0445023148148148148,
    -0.0259847148736037625,
    +0.0156356325323339212,
    -0.00961689202429943171,
    +0.00601454325295611786,
    -0.00381129803489199923,
    +0.00244087799114398267,
    -0.00157693034468678425,
    +0.00102626332050760715,
    -0.000672061631156136204,
    +0.000442473061814620910,
    -0.000292677224729627445,
    +0.000194387276054539318,
    -0.000129574266852748819,
    +0.0000866503580520812717,
    -0.0000581136075044138168
  };
  const double ap = abs(p);
  if (ap < 0.01159)
    return
      -1 +
      p*(1 +
      p*(q[2] +
      p*(q[3] +
      p*(q[4] +
      p*(q[5] +
      p*q[6]
      )))));
  else if (ap < 0.0766)
    return
      -1 +
      p*(1 +
      p*(q[2] +
      p*(q[3] +
      p*(q[4] +
      p*(q[5] +
      p*(q[6] +
      p*(q[7] +
      p*(q[8] +
      p*(q[9] +
      p*q[10]
      )))))))));
  else
    return
      -1 +
      p*(1 +
      p*(q[2] +
      p*(q[3] +
      p*(q[4] +
      p*(q[5] +
      p*(q[6] +
      p*(q[7] +
      p*(q[8] +
      p*(q[9] +
      p*(q[10] +
      p*(q[11] +
      p*(q[12] +
      p*(q[13] +
      p*(q[14] +
      p*(q[15] +
      p*(q[16] +
      p*(q[17] +
      p*(q[18] +
      p*(q[19] +
      p*q[20]
      )))))))))))))))))));
}


inline double LambertW0ZeroSeries(const double z)
{
  return
    z*(1 -
    z*(1 -
    z*(1.5 -
    z*(2.6666666666666666667 -
    z*(5.2083333333333333333 -
    z*(10.8 -
    z*(23.343055555555555556 -
    z*(52.012698412698412698 -
    z*(118.62522321428571429 -
    z*(275.57319223985890653 -
    z*(649.78717234347442681 -
    z*(1551.1605194805194805 -
    z*(3741.4497029592385495 -
    z*(9104.5002411580189358 -
    z*(22324.308512706601434 -
    z*(55103.621972903835338 -
    z*136808.86090394293563
    ))))))))))))))));
}


inline double FinalResult(const double w, const double y)
{
  const double f0 = w - y;
  const double f1 = 1 + y;
  const double f00 = f0 * f0;
  const double f11 = f1 * f1;
  const double f0y = f0 * y;
  return w - 4 * f0 * (6 * f1 * (f11 + f0y) + f00 * y) /
                      (f11 * (24 * f11 + 36 * f0y) +
                       f00 * (6 * y * y + 8 * f1 * y + f0y));
}


double Randomizer::LambertW0(const double z)
{
  static double e[66];
  static double g[65];
  static double a[12];
  static double b[12];

  if (!e[0]) {
    const double e1 = 1 / M_E;
    double ej = 1;
    e[0] = M_E;
    e[1] = 1;
    g[0] = 0;
    for (int j = 1, jj = 2; jj < 66; ++jj) {
      ej *= M_E;
      e[jj] = e[j] * e1;
      g[j] = j * ej;
      j = jj;
    }
    a[0] = std::sqrt(e1);
    b[0] = 0.5;
    for (int j = 0, jj = 1; jj < 12; ++jj) {
      a[jj] = std::sqrt(a[j]);
      b[jj] = b[j] * 0.5;
      j = jj;
    }
  }
  if (std::abs(z) < 0.05)
    return LambertW0ZeroSeries(z);
  if (z < -0.35) {
    const double p2 = 2 * (M_E * z + 1);
    if (p2 > 0)
      return LambertWSeries(std::sqrt(p2));
    if (p2 == 0)
      return -1;
    throw std::runtime_error("(lambertw0) Argument out of range.");
  }
  int n;
  for (n = 0; n <= 2; ++n)
    if (g[n] > z)
      goto line1;
  n = 2;
  for (int j = 1; j <= 5; ++j) {
    n *= 2;
    if (g[n] > z)
      goto line2;
  }
  throw std::runtime_error("(lambertw0) Argument too large.");
line2:
  {
    int nh = n / 2;
    for (int j = 1; j <= 5; ++j) {
      nh /= 2;
      if (nh <= 0)
        break;
      if (g[n-nh] > z)
        n -= nh;
    }
  }
line1:
  --n;
  int jmax = 8;
  if (z <= -0.36)
    jmax = 12;
  else if (z <= -0.3)
    jmax = 11;
  else if (n <= 0)
    jmax = 10;
  else if (n <= 1)
    jmax = 9;
  double y = z * e[n+1];
  double w = n;
  for (int j = 0; j < jmax; ++j) {
    const double wj = w + b[j];
    const double yj = y * a[j];
    if (wj < yj) {
      w = wj;
      y = yj;
    }
  }
  return FinalResult(w, y);
}
// End of Lambert W function implementation

// Values of log(lambda) to use for foundress-Poisson distribution, starting from a = 1/256 to a = 255/256 in steps of 1/256.
double Randomizer::FoundressPoisson_LogLambda[256] = {
    5.5490914045946837518, 4.8598739376714563676, 4.4583548216393218411, 4.1746355810215849402, 3.9554723034763945577,
    3.7771491107553165634, 3.6270155863428330534, 3.4975209021982966995, 3.3837949606719561757, 3.2825128415521391823,
    3.1913033792074050332, 3.1084161826257843408, 3.0325224056294168840, 2.9625895584403929561, 2.8977995927504300866,
    2.8374934331027250600, 2.7811322902699306958, 2.7282699772456133758, 2.6785326529953934482, 2.6316037173053619114,
    2.5872123676509315437, 2.5451248185422103987, 2.5051374946038991176, 2.4670717101882773115, 2.4307694821294680843,
    2.3960902135522448297, 2.3629080508791155957, 2.3311097627848256231, 2.3005930246189114641, 2.2712650183267588666,
    2.2430412783953856959, 2.2158447303007200446, 2.1896048803531127369, 2.1642571254722251517, 2.1397421588484801802,
    2.1160054531285581447, 2.0929968070738049768, 2.0706699448896435101, 2.0489821598590083340, 2.0278939957305648356,
    2.0073689606678595254, 1.9873732695777535096, 1.9678756113962181384, 1.9488469384847255661, 1.9302602757311730919,
    1.9120905472922897772, 1.8943144191862868464, 1.8769101561641798881, 1.8598574914693211113, 1.8431375082468526294,
    1.8267325314949027781, 1.8106260295621630085, 1.7948025242955778502, 1.7792475090296258067, 1.7639473736870272536,
    1.7488893363309332418, 1.7340613805718745333, 1.7194521982897528201, 1.7050511371826571061, 1.6908481527008452083,
    1.6768337639662702632, 1.6629990133161332011, 1.6493354291433091063, 1.6358349917375756277, 1.6224901018596895597,
    1.6092935518056998845, 1.5962384987418021254, 1.5833184401107398553, 1.5705271909293907484, 1.5578588628140852546,
    1.5453078445853500877, 1.5328687843175468064, 1.5205365727112536423, 1.5083063276774182349, 1.4961733800324275023,
    1.4841332602123664675, 1.4721816859229783780, 1.4603145506492989369, 1.4485279129556871691, 1.4368179865130528139,
    1.4251811307956376851, 1.4136138423946635889, 1.4021127469007348409, 1.3906745913109561652, 1.3792962369204679884,
    1.3679746526614877666, 1.3567069088560006485, 1.3454901713510472039, 1.3343216960080948041, 1.3231988235202678528,
    1.3121189745333428078, 1.3010796450483075315, 1.2900784020850375455, 1.2791128795882311664, 1.2681807745581985536,
    1.2572798433904310844, 1.2464078984090793956, 1.2355628045805810977, 1.2247424763946948012, 1.2139448749011210182,
    1.2031680048907353875, 1.1924099122112450377, 1.1816686812077812352, 1.1709424322796015971, 1.1602293195446662377,
    1.1495275286044117635, 1.1388352744015306506, 1.1281507991640415955, 1.1174723704293623161, 1.1067982791424477718,
    1.0961268378224759967, 1.0854563787928213081, 1.0747852524693937148, 1.0641118257026898064, 1.0534344801691450932,
    1.0427516108075869372, 1.0320616242968407850, 1.0213629375706771985, 1.0106539763664943443, 0.99993317380427182428,
    0.98919896899246673172, 0.97844980565765471425, 0.96768413079483150963, 0.95690039333536858202, 0.94609704282972795220,
    0.93527252814210193321, 0.92442529615420121480, 0.91355379047546680926, 0.90265645015704287779, 0.89173170840684001792,
    0.88077799130307732334, 0.86979371650367520719, 0.85877729194888252628, 0.84772711455451466023, 0.83664156889315188792,
    0.82551902586063019562, 0.81435784132510335986, 0.80315635475592217496, 0.79191288782950652880, 0.78062574300931997229,
    0.76929320209697371613, 0.75791352475140416622, 0.74648494697294431877, 0.73500567954900697387, 0.72347390645798004272,
    0.71188778322774348695, 0.70024543524511828618, 0.68854495601233656199, 0.67678440534645978310, 0.66496180751742528514,
    0.65307514932023225107, 0.64112237807646343946, 0.62910139956010391327, 0.61701007584231570835, 0.60484622304950252936,
    0.59260760902863973687, 0.58029195091348084734, 0.56789691258482033476, 0.55542010201753710952, 0.54285906850668952384,
    0.53021129976434855369, 0.51747421887829947451, 0.50464518112312062303, 0.49172147061344534391, 0.47870029678844322474,
    0.46557879071579494346, 0.45235400120247148958, 0.43902289069870409355, 0.42558233098045000764, 0.41202909859448888508,
    0.39835987004899742203, 0.38457121673108879412, 0.37065959953123667203, 0.35662136315282488841, 0.34245273008322496544,
    0.32814979420076678673, 0.31370851398972732227, 0.29912470533299456710, 0.28439403384931910557, 0.26951200673909808669,
    0.25447396409927247607, 0.23927506966423986445, 0.22391030092558442122, 0.20837443857889523247, 0.19266205524092211432,
    0.17676750337459781748, 0.16068490235329987992, 0.14440812458862692380, 0.12793078063816637480, 0.11124620320095827963,
    0.094347429898460274944, 0.077227184727783088070, 0.059877858061459093841, 0.042291485053985305997, 0.024459722299361594255, 
    0.0063738225660197678560, -0.011975392585108220189, -0.030597562516649728231, -0.049502819808872750018, -0.068701825405796557167, 
    -0.088205806921918511465, -0.10802660045778708642, -0.12817669631822678489, -0.14866928907951046379, -0.16951833251242237610, 
    -0.19073859993853778394, -0.21234575067854272179, -0.23435640334665591711, -0.25678821685632552407, -0.27965998013256737620,
    -0.30299171167965555096, -0.32680477033341337467, -0.35112197874149247978, -0.37596776136922083200, -0.40136829913192989538, 
    -0.42735170311793696518, -0.45394821030318399657, -0.48119040468621782081, -0.50911346791255374100, -0.53775546423799680529, 
    -0.56715766563744574036, -0.59736492404426044800, -0.62842609916568081818, -0.66039455213985598370, -0.69332871758270142593, 
    -0.72729276945440679558, -0.76235739983949002418, -0.79860073442730838966, -0.83610941453945741841, -0.87497988344137489491,
    -0.91531992504392012400, -0.95725051685849371630, -1.0009080775167347177, -1.0464472141795768678, -1.0940441094954531653, 
    -1.1439007355214256823, -1.1962501493984540879, -1.2513632221395034616, -1.3095572926906466904, -1.3712074486477998647, 
    -1.4367614524936178633, -1.5067598253616307780, -1.5818633866598719173, -1.6628918393348275373, -1.7508791842965669705, 
    -1.8471556211944450965, -1.9534727482714073776, -2.0722028139462080887, -2.2066717185478532670, -2.3617504547707843798, 
    -2.5449906739549490453, -2.7690435543599964952, -3.0576256824558680769, -3.4639816317729938966, -4.1580104972824383225, -1000
};
