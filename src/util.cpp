#include "util.h"

namespace Util {

    std::vector<double> generateBinomialCoefficients(int n) {
        std::vector<double> coeff(n+1);

        for (int i = 0; i < n+1; i++) {
            coeff.push_back(C(n, i));
        }
        
        return coeff;
    }

    /**
     * https://cp-algorithms.com/combinatorics/binomial-coefficients.html#improved-implementation
     * //TODO could cache for repeated computations.. always computing (n i) for i in [0, n]
    */
    int C(int n, int k) {
        double res = 1;
        for (int i = 1; i <= k; ++i)
            res = res * (n - k + i) / i;
        return (int)(res + 0.01);
    }
}