#pragma once

#include <vector>

namespace Util {

    /**
     * @brief Generates binomial coefficients [(n 0), (n 1), ..., (n n)]
    */
    std::vector<double> generateBinomialCoefficients(int n);

}