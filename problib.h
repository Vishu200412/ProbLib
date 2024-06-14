#ifndef PROBLIB_H
#define PROBLIB_H

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <random>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <chrono>

namespace ProbLib {
    // Geometric Distribution
    class GeometricDistribution {
    public:
        explicit GeometricDistribution(double probability)
            : p(probability) {
            if (p <= 0.0 || p >= 1.0) {
                throw std::invalid_argument("Probability must be between 0 and 1 (exclusive).");
            }
        }

        double pmf(int k) const {
            if (k <= 0) {
                throw std::invalid_argument("Iteration (k) must be a positive integer.");
            }
            return p * std::pow(1 - p, k - 1);
        }

        double cdf(int k) const {
            if (k <= 0) {
                throw std::invalid_argument("Iteration (k) must be a positive integer.");
            }
            return 1 - std::pow(1 - p, k);
        }

        double mean() const {
            return 1 / p;
        }

        double variance() const {
            return (1 - p) / (p * p);
        }

    private:
        double p;
    };

    // Poisson Distribution
    class PoissonDistribution {
    public:
        explicit PoissonDistribution(double lambda)
            : lambda(lambda) {
            if (lambda < 0.0) {
                throw std::invalid_argument("Lambda must be non-negative.");
            }
        }

        double pmf(int k) const {
            if (k < 0) {
                throw std::invalid_argument("Iteration (k) must be a non-negative integer.");
            }
            return std::exp(-lambda) * std::pow(lambda, k) / std::tgamma(k + 1);
        }

        double cdf(int k) const {
            if (k < 0) {
                throw std::invalid_argument("Iteration (k) must be a non-negative integer.");
            }
            double sum = 0.0;
            for (int i = 0; i <= k; ++i) {
                sum += pmf(i);
            }
            return sum;
        }

        double mean() const {
            return lambda;
        }

        double variance() const {
            return lambda;
        }

    private:
        double lambda;
    };

    // Binomial Distribution
    class BinomialDistribution {
    public:
        BinomialDistribution(int trials, double probability)
            : n(trials), p(probability) {
            if (n <= 0) {
                throw std::invalid_argument("Number of trials (n) must be a positive integer.");
            }
            if (p <= 0.0 || p >= 1.0) {
                throw std::invalid_argument("Probability must be between 0 and 1 (exclusive).");
            }
        }

        double pmf(int k) const {
            if (k < 0 || k > n) {
                throw std::invalid_argument("Invalid value for k.");
            }
            double coeff = std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
            return coeff * std::pow(p, k) * std::pow(1 - p, n - k);
        }

        double cdf(int k) const {
            if (k < 0 || k > n) {
                throw std::invalid_argument("Invalid value for k.");
            }
            double sum = 0.0;
            for (int i = 0; i <= k; ++i) {
                sum += pmf(i);
            }
            return sum;
        }

        double mean() const {
            return n * p;
        }

        double variance() const {
            return n * p * (1 - p);
        }

    private:
        int n;
        double p;
    };

    // Exponential Distribution
    class ExponentialDistribution {
    public:
        explicit ExponentialDistribution(double mean)
            : mean(mean), lambda(1.0 / mean) {
            if (mean <= 0) {
                throw std::invalid_argument("Mean must be positive.");
            }
        }

        double pdf(double x) const {
            if (x < 0) {
                throw std::invalid_argument("x must be non-negative.");
            }
            return lambda * std::exp(-lambda * x);
        }

        double cdf(double x) const {
            if (x < 0) {
                throw std::invalid_argument("x must be non-negative.");
            }
            return 1.0 - std::exp(-lambda * x);
        }

        double prob_gte(double a) const {
            if (a < 0) {
                throw std::invalid_argument("a must be non-negative.");
            }
            return std::exp(-lambda * a);
        }

        double prob_range(double a, double b) const {
            if (a < 0 || b < 0) {
                throw std::invalid_argument("a and b must be non-negative.");
            }
            if (b < a) {
                throw std::invalid_argument("b must be greater than or equal to a.");
            }
            return cdf(b) - cdf(a);
        }

        double get_mean() const {
            return mean;
        }

        double variance() const {
            return 1.0 / (lambda * lambda);
        }

    private:
        double mean;
        double lambda;
    };

    // Normal Distribution
    class NormalDistribution {
    public:
        NormalDistribution(double mean, double variance)
            : mean(mean), variance(variance), stddev(std::sqrt(variance)) {
            if (variance <= 0) {
                throw std::invalid_argument("Variance must be positive.");
            }
        }

        double pdf(double x) const {
            return (1.0 / (stddev * std::sqrt(2 * M_PI))) * std::exp(-0.5 * std::pow((x - mean) / stddev, 2));
        }

        double cdf(double x) const {
            return 0.5 * (1 + std::erf((x - mean) / (stddev * std::sqrt(2))));
        }

        double prob_range(double a, double b) const {
            if (b < a) {
                throw std::invalid_argument("b must be greater than or equal to a.");
            }
            return cdf(b) - cdf(a);
        }

        double get_mean() const {
            return mean;
        }

        double get_variance() const {
            return variance;
        }

    private:
        double mean;
        double variance;
        double stddev;
    };

    // Markov Chain class for discrete Markov chains
    class MarkovChain {
    public:
        explicit MarkovChain(const std::vector<std::vector<double>>& transition_matrix)
            : transition_matrix(transition_matrix) {
            validate_transition_matrix();
        }
        
        // Helper function to check if the Markov chain is irreducible
        bool is_irreducible(const std::vector<std::vector<double>>& transition_matrix) {
            size_t n = transition_matrix.size();
            std::vector<std::unordered_set<int>> reachability(n);
            for (size_t i = 0; i < n; ++i) {
                std::queue<int> q;
                q.push(i);
                while (!q.empty()) {
                    int current = q.front();
                    q.pop();
                    for (size_t j = 0; j < n; ++j) {
                        if (transition_matrix[current][j] > 0 && reachability[i].find(j) == reachability[i].end()) {
                            reachability[i].insert(j);
                            q.push(j);
                        }
                    }
                }
            }
            for (size_t i = 0; i < n; ++i) {
                if (reachability[i].size() != n) {
                    return false;
                }
            }
            return true;
        }


        // Helper function to check if the Markov chain is aperiodic
        bool is_aperiodic(const std::vector<std::vector<double>>& transition_matrix) {
        size_t n = transition_matrix.size();
        std::vector<int> periods(n, -1);
        std::vector<bool> visited(n, false);

        for (size_t i = 0; i < n; ++i) {
            if (periods[i] == -1) {
            periods[i] = 0;
            std::queue<int> q;
            q.push(i);
            visited[i] = true;

            while (!q.empty()) {
                int current = q.front();
                q.pop();

                for (size_t j = 0; j < n; ++j) {
                if (transition_matrix[current][j] > 0 && !visited[j]) {
                    if (periods[j] == -1) {
                    periods[j] = periods[current] + 1;
                    q.push(j);
                    visited[j] = true;
                    } else {
                    if ((periods[j] - periods[current]) % 2 != 0) {
                        return true;
                    }
                    }
                }
                }
            }
            }
        }
        return false;
        }


        // Function to check if the Markov chain is ergodic
        bool is_ergodic(const std::vector<std::vector<double>>& transition_matrix) {
            return is_irreducible(transition_matrix) && is_aperiodic(transition_matrix);
        }

        // Calculate n-step transition probabilities using Chapman-Kolmogorov equation
        std::vector<std::vector<double>> n_step_transition(int n) const {
            if (n < 1) {
                throw std::invalid_argument("Number of steps must be a positive integer.");
            }

            std::vector<std::vector<double>> result = transition_matrix;
            for (int step = 1; step < n; ++step) {
                result = matrix_multiply(result, transition_matrix);
            }
            return result;
        }

        // Function to compute the steady-state distribution
        std::vector<std::vector<double>> compute_steady_state_matrix(const std::vector<std::vector<double>>& transition_matrix, int max_iterations = 1000, double tolerance = 1e-10, double time_limit_seconds = 5.0) {
            size_t n = transition_matrix.size();
            std::vector<std::vector<double>> steady_state_matrix(n, std::vector<double>(n, 0.0));

            auto start_time = std::chrono::high_resolution_clock::now();

            for (size_t initial_state = 0; initial_state < n; ++initial_state) {
                std::vector<double> distribution(n, 0.0);
                distribution[initial_state] = 1.0;

                for (int iter = 0; iter < max_iterations; ++iter) {
                    std::vector<double> next_distribution(n, 0.0);
                    for (size_t i = 0; i < n; ++i) {
                        for (size_t j = 0; j < n; ++j) {
                            next_distribution[i] += distribution[j] * transition_matrix[j][i];
                        }
                    }

                    // Check for convergence
                    double diff = 0.0;
                    for (size_t i = 0; i < n; ++i) {
                        diff += std::abs(next_distribution[i] - distribution[i]);
                    }
                    if (diff < tolerance) {
                        break;
                    }

                    distribution = next_distribution;

                    // Check for time limit
                    auto current_time = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = current_time - start_time;
                    if (elapsed.count() > time_limit_seconds) {
                        throw std::runtime_error("Matrix is too large to compute the steady-state distribution within the time limit.");
                    }
                }

                steady_state_matrix[initial_state] = distribution;
            }

            return steady_state_matrix;
        }

        // Wrapper function to compute the steady-state distribution or steady-state matrix
        std::pair<std::string, std::vector<std::vector<double>>> steady_state(const std::vector<std::vector<double>>& transition_matrix) {
            if(!is_aperiodic(transition_matrix)){
                return {"Markov chaini is periodic, probability keeps oscillating",{{}}};
            }
            else if (is_ergodic(transition_matrix)) {
                // Ergodic case: compute a single steady-state distribution
                std::vector<double> distribution(transition_matrix.size(), 1.0 / transition_matrix.size());
                std::vector<double> next_distribution = distribution;
                const double tolerance = 1e-10;
                do {
                    distribution = next_distribution;
                    for (size_t i = 0; i < transition_matrix.size(); ++i) {
                        next_distribution[i] = 0.0;
                        for (size_t j = 0; j < transition_matrix.size(); ++j) {
                            next_distribution[i] += distribution[j] * transition_matrix[j][i];
                        }
                    }
                } while (std::inner_product(distribution.begin(), distribution.end(), next_distribution.begin(), 0.0, std::plus<double>(), [](double a, double b) { return std::abs(a - b); }) > tolerance);

                std::vector<std::vector<double>> result_matrix(1, next_distribution);
                return { "Markov chain is ergodic.", result_matrix };
            }
            else {
                // Non-ergodic and aperiodic case: compute the steady-state distribution matrix
                return { "Markov chain is not ergodic. Steady-state distribution is not unique.", compute_steady_state_matrix(transition_matrix) };
            }
        }


    private:
        std::vector<std::vector<double>> transition_matrix;

        void validate_transition_matrix() const {
            for (const auto& row : transition_matrix) {
                if (row.size() != transition_matrix.size()) {
                    throw std::invalid_argument("Transition matrix must be square.");
                }
                double sum = std::accumulate(row.begin(), row.end(), 0.0);
                if (std::abs(sum - 1.0) > 1e-10) {
                    throw std::invalid_argument("Each row of the transition matrix must sum to 1.");
                }
            }
        }

        // Helper function to multiply two matrices
        std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) const {
            size_t n = A.size();
            std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    for (size_t k = 0; k < n; ++k) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            return C;
        }
    };
}

#endif