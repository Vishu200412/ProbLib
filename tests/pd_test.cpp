#include <iostream>
#include <cassert>
#include "problib.h"

void testGeo() {
    try {
        ProbLib::GeometricDistribution geoDist(0.5);
        double value = geoDist.pmf(3);
        assert(value == 0.125);
        std::cout << "Geo PMF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Geo PMF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::GeometricDistribution geoDist(1.5); // Should throw an exception
        std::cerr << "Geo PMF test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Geo PMF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::GeometricDistribution geoDist(0.5);
        geoDist.pmf(-3); // Should throw an exception
        std::cerr << "Geo PMF test failed: No exception thrown for invalid k" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Geo PMF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testGeomean() {
    try {
        ProbLib::GeometricDistribution geoDist(0.5);
        double value = geoDist.mean();
        assert(value == 2);
        std::cout << "Geo Mean test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Geo Mean test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::GeometricDistribution geoDist(-0.5); // Should throw an exception
        std::cerr << "Geo Mean test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Geo Mean test passed with expected exception: " << e.what() << std::endl;
    }
}

void testGeovar() {
    try {
        ProbLib::GeometricDistribution geoDist(0.5);
        double value = geoDist.variance();
        assert(value == 2);
        std::cout << "Geo Variance test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Geo Variance test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::GeometricDistribution geoDist(1.5); // Should throw an exception
        std::cerr << "Geo Variance test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Geo Variance test passed with expected exception: " << e.what() << std::endl;
    }
}

void testGeocdf() {
    try {
        ProbLib::GeometricDistribution geoDist(0.5);
        double value = geoDist.cdf(3);
        assert(value == 0.875);
        std::cout << "Geo CDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Geo CDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::GeometricDistribution geoDist(1.5); // Should throw an exception
        std::cerr << "Geo CDF test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Geo CDF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::GeometricDistribution geoDist(0.5);
        geoDist.cdf(-3); // Should throw an exception
        std::cerr << "Geo CDF test failed: No exception thrown for invalid n" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Geo CDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testPoisson() {
    try {
        ProbLib::PoissonDistribution poissonDist(4.0);
        double value = poissonDist.pmf(3);
        assert(value > 0);
        std::cout << "Poisson PMF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Poisson PMF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::PoissonDistribution poissonDist(-4.0); // Should throw an exception
        std::cerr << "Poisson PMF test failed: No exception thrown for invalid lambda" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Poisson PMF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::PoissonDistribution poissonDist(4.0);
        poissonDist.pmf(-3); // Should throw an exception
        std::cerr << "Poisson PMF test failed: No exception thrown for invalid k" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Poisson PMF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testPoissonVar() {
    try {
        ProbLib::PoissonDistribution poissonDist(4.0);
        double value = poissonDist.variance();
        assert(value == 4.0);
        std::cout << "Poisson Variance test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Poisson Variance test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::PoissonDistribution poissonDist(-4.0); // Should throw an exception
        std::cerr << "Poisson Variance test failed: No exception thrown for invalid lambda" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Poisson Variance test passed with expected exception: " << e.what() << std::endl;
    }
}

void testPoissonCdf() {
    try {
        ProbLib::PoissonDistribution poissonDist(4.0);
        double value = poissonDist.cdf(3);
        assert(value >= 0 && value <= 1);
        std::cout << "Poisson CDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Poisson CDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::PoissonDistribution poissonDist(-4.0); // Should throw an exception
        std::cerr << "Poisson CDF test failed: No exception thrown for invalid lambda" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Poisson CDF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::PoissonDistribution poissonDist(4.0);
        poissonDist.cdf(-3); // Should throw an exception
        std::cerr << "Poisson CDF test failed: No exception thrown for invalid k" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Poisson CDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testBinomial() {
    try {
        ProbLib::BinomialDistribution binomDist(4, 0.5);
        double value = binomDist.pmf(3);
        assert(value == 0.25);
        std::cout << "Binomial PMF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Binomial PMF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(10, 0.5);
        binomDist.pmf(-3); // Should throw an exception
        std::cerr << "Binomial PMF test failed: No exception thrown for invalid k" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial PMF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(-10, 0.5); // Should throw an exception
        std::cerr << "Binomial PMF test failed: No exception thrown for invalid n" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial PMF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(10, 1.5); // Should throw an exception
        std::cerr << "Binomial PMF test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial PMF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testBinomialMean() {
    try {
        ProbLib::BinomialDistribution binomDist(10, 0.5);
        double value = binomDist.mean();
        assert(value == 5);
        std::cout << "Binomial Mean test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Binomial Mean test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(-10, 0.5); // Should throw an exception
        std::cerr << "Binomial Mean test failed: No exception thrown for invalid n" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial Mean test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(10, 1.5); // Should throw an exception
        std::cerr << "Binomial Mean test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial Mean test passed with expected exception: " << e.what() << std::endl;
    }
}

void testBinomialVar() {
    try {
        ProbLib::BinomialDistribution binomDist(10, 0.5);
        double value = binomDist.variance();
        assert(value == 2.5);
        std::cout << "Binomial Variance test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Binomial Variance test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(-10, 0.5); // Should throw an exception
        std::cerr << "Binomial Variance test failed: No exception thrown for invalid n" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial Variance test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(10, 1.5); // Should throw an exception
        std::cerr << "Binomial Variance test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial Variance test passed with expected exception: " << e.what() << std::endl;
    }
}

void testBinomialCdf() {
    try {
        ProbLib::BinomialDistribution binomDist(10, 0.5);
        double value = binomDist.cdf(3);
        assert(value >= 0 && value <= 1);
        std::cout << "Binomial CDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Binomial CDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(-10, 0.5); // Should throw an exception
        std::cerr << "Binomial CDF test failed: No exception thrown for invalid n" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial CDF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(10, 1.5); // Should throw an exception
        std::cerr << "Binomial CDF test failed: No exception thrown for invalid probability" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial CDF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::BinomialDistribution binomDist(10, 0.5);
        binomDist.cdf(-3); // Should throw an exception
        std::cerr << "Binomial CDF test failed: No exception thrown for invalid k" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Binomial CDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testExponentialPdf() {
    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        double value = expDist.pdf(1.0);
        assert(value >= 0);
        std::cout << "Exponential PDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exponential PDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(-2.0); // Should throw an exception
        std::cerr << "Exponential PDF test failed: No exception thrown for invalid mean" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential PDF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        expDist.pdf(-1.0); // Should throw an exception
        std::cerr << "Exponential PDF test failed: No exception thrown for invalid x" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential PDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testExponentialMean() {
    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        double value = expDist.get_mean();
        assert(value == 2.0);
        std::cout << "Exponential Mean test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exponential Mean test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(-2.0); // Should throw an exception
        std::cerr << "Exponential Mean test failed: No exception thrown for invalid mean" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Mean test passed with expected exception: " << e.what() << std::endl;
    }
}

void testExponentialVar() {
    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        double value = expDist.variance();
        assert(value >= 0);
        std::cout << "Exponential Variance test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exponential Variance test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(-2.0); // Should throw an exception
        std::cerr << "Exponential Variance test failed: No exception thrown for invalid mean" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Variance test passed with expected exception: " << e.what() << std::endl;
    }
}

void testExponentialProbGte() {
    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        double value = expDist.prob_gte(1.0);
        assert(value >= 0 && value <= 1);
        std::cout << "Exponential Probability GTE test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exponential Probability GTE test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(-2.0); // Should throw an exception
        std::cerr << "Exponential Probability GTE test failed: No exception thrown for invalid mean" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Probability GTE test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        expDist.prob_gte(-1.0); // Should throw an exception
        std::cerr << "Exponential Probability GTE test failed: No exception thrown for invalid a" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Probability GTE test passed with expected exception: " << e.what() << std::endl;
    }
}

void testExponentialCdf() {
    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        double value = expDist.cdf(1.0);
        assert(value >= 0 && value <= 1);
        std::cout << "Exponential CDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exponential CDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(-2.0); // Should throw an exception
        std::cerr << "Exponential CDF test failed: No exception thrown for invalid mean" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential CDF test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        expDist.cdf(-1.0); // Should throw an exception
        std::cerr << "Exponential CDF test failed: No exception thrown for invalid a" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential CDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testExpProbRange() {
    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        double value = expDist.prob_range(1.0, 3.0);
        assert(value >= 0 && value <= 1);
        std::cout << "Exponential Prob Range test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exponential Prob Range test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(-2.0); // Should throw an exception
        std::cerr << "Exponential Prob Range test failed: No exception thrown for invalid mean" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Prob Range test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        expDist.prob_range(-1.0, 3.0); // Should throw an exception
        std::cerr << "Exponential Prob Range test failed: No exception thrown for invalid a" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Prob Range test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        expDist.prob_range(1.0, -3.0); // Should throw an exception
        std::cerr << "Exponential Prob Range test failed: No exception thrown for invalid b" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Prob Range test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::ExponentialDistribution expDist(2.0);
        expDist.prob_range(3.0, 1.0); // Should throw an exception
        std::cerr << "Exponential Prob Range test failed: No exception thrown for b < a" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Exponential Prob Range test passed with expected exception: " << e.what() << std::endl;
    }
}

void testNormalPdf() {
    try {
        ProbLib::NormalDistribution normDist(0.0, 1.0);
        double value = normDist.pdf(1.0);
        assert(value > 0);
        std::cout << "Normal PDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Normal PDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::NormalDistribution normDist(0.0, -1.0); // Should throw an exception
        std::cerr << "Normal PDF test failed: No exception thrown for invalid variance" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Normal PDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testNormalCdf() {
    try {
        ProbLib::NormalDistribution normDist(0.0, 1.0);
        double value = normDist.cdf(1.0);
        assert(value >= 0 && value <= 1);
        std::cout << "Normal CDF test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Normal CDF test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::NormalDistribution normDist(0.0, -1.0); // Should throw an exception
        std::cerr << "Normal CDF test failed: No exception thrown for invalid variance" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Normal CDF test passed with expected exception: " << e.what() << std::endl;
    }
}

void testNormalProbRange() {
    try {
        ProbLib::NormalDistribution normDist(0.0, 1.0);
        double value = normDist.prob_range(-1.0, 1.0);
        assert(value >= 0 && value <= 1);
        std::cout << "Normal Prob Range test passed with value: " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Normal Prob Range test failed: " << e.what() << std::endl;
    }

    try {
        ProbLib::NormalDistribution normDist(0.0, -1.0); // Should throw an exception
        std::cerr << "Normal Prob Range test failed: No exception thrown for invalid variance" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Normal Prob Range test passed with expected exception: " << e.what() << std::endl;
    }

    try {
        ProbLib::NormalDistribution normDist(0.0, 1.0);
        normDist.prob_range(1.0, -1.0); // Should throw an exception
        std::cerr << "Normal Prob Range test failed: No exception thrown for invalid range" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Normal Prob Range test passed with expected exception: " << e.what() << std::endl;
    }
}

int main() {
    testGeo();
    testGeomean();
    testGeovar();
    testGeocdf();

    testPoisson();
    testPoissonVar();
    testPoissonCdf();

    testBinomial();
    testBinomialMean();
    testBinomialVar();
    testBinomialCdf();

    testExponentialPdf();
    testExponentialMean();
    testExponentialVar();
    testExponentialProbGte();
    testExponentialCdf();
    testExpProbRange();

    testNormalPdf();
    testNormalCdf();
    testNormalProbRange();

    return 0;
}
