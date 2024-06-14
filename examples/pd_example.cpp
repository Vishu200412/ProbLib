#include <iostream>
#include <exception>
#include "problib.h"

int main() {
    try {
        double p = 0.5;
        int k = 3;
        int n = 10;

        // Geometric distribution
        ProbLib::GeometricDistribution geoDist(p);

        // Input probability of success (p) in one trial and number of trials to get first success (k) 
        // to get probability that first success occurs at kth trial
        double geo_pmf = geoDist.pmf(k);
        std::cout << "Geo PMF: " << geo_pmf << std::endl;

        // Mean number of trials to get success
        double geo_mean = geoDist.mean();
        std::cout << "Geo Mean: " << geo_mean << std::endl;

        // Variance
        double geo_var = geoDist.variance();
        std::cout << "Geo Variance: " << geo_var << std::endl;

        // P(X ≤ k); probability success occurs within k trials
        double geo_cdf = geoDist.cdf(k);
        std::cout << "Geo CDF: " << geo_cdf << std::endl;

        // Poisson distribution
        // We are always given in question the mean (lambda parameter)
        double mean = 4.0;
        ProbLib::PoissonDistribution poissonDist(mean);

        // PMF of Poisson distribution with mean value of no. of occurrences per unit time as 'mean' and no. of occurrences k in unit time
        double poisson_pmf = poissonDist.pmf(k);
        std::cout << "Poisson PMF: " << poisson_pmf << std::endl;

        // Variance
        double poisson_var = poissonDist.variance();
        std::cout << "Poisson Variance: " << poisson_var << std::endl;

        // P(X ≤ k); (CDF) of the Poisson distribution represents the probability that the number of occurrences in a given interval is less than or equal to k
        double poisson_cdf = poissonDist.cdf(k);
        std::cout << "Poisson CDF: " << poisson_cdf << std::endl;

        // Binomial distribution
        int binom_n = 10;
        double binom_p = 0.3;
        ProbLib::BinomialDistribution binomDist(binom_n, binom_p);

        // Probability there will be k successful outcomes out of binom_n trials with success in each trial having probability of binom_p;
        double binom_pmf = binomDist.pmf(k);
        std::cout << "Binomial PMF: " << binom_pmf << std::endl;

        // Mean no. of successful outcomes in n trials
        double binom_mean = binomDist.mean();
        std::cout << "Binomial Mean: " << binom_mean << std::endl;

        // Variance
        double binom_var = binomDist.variance();
        std::cout << "Binomial Variance: " << binom_var << std::endl;

        // CDF; P(X ≤ k); probability that no. of successful outcomes is less than or equal to k
        double binom_cdf = binomDist.cdf(k);
        std::cout << "Binomial CDF: " << binom_cdf << std::endl;

        // Exponential distribution
        mean = 2.0;
        double x = 1.0;
        ProbLib::ExponentialDistribution expDist(mean);

        // PDF
        double exp_pdf = expDist.pdf(x);
        std::cout << "Exponential PDF: " << exp_pdf << std::endl;

        // Mean
        double exp_mean = expDist.get_mean();
        std::cout << "Exponential Mean: " << exp_mean << std::endl;

        // Variance
        double exp_var = expDist.variance();
        std::cout << "Exponential Variance: " << exp_var << std::endl;

        // P(X >= a)
        double exp_prob_gte = expDist.prob_gte(x);
        std::cout << "Exponential Prob(x >= a): " << exp_prob_gte << std::endl;

        // CDF
        double exp_cdf = expDist.cdf(x);
        std::cout << "Exponential CDF: " << exp_cdf << std::endl;

        // Probability in range [a, b]
        double exp_prob_range = expDist.prob_range(1.0, 3.0);
        std::cout << "Exponential Prob(1 <= x <= 3): " << exp_prob_range << std::endl;

        // Normal distribution
        double normal_mean = 0.0;
        double normal_variance = 1.0;
        ProbLib::NormalDistribution normDist(normal_mean, normal_variance);

        // PDF
        double normal_pdf = normDist.pdf(x);
        std::cout << "Normal PDF: " << normal_pdf << std::endl;

        // Probability that a normally distributed random variable X takes on a value less than or equal to a specific value x.
        double normal_cdf = normDist.cdf(x);
        std::cout << "Normal CDF: " << normal_cdf << std::endl;

        // Calculating the probability that x takes a value between a specified range
        double normal_prob_range = normDist.prob_range(-1.0, 1.0);
        std::cout << "Normal Prob(-1 <= x <= 1): " << normal_prob_range << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
