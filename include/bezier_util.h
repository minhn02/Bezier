#pragma once

#include <vector>
#include <map>

namespace Util {

    /**
     * @brief Returns binomial coefficients [(n 0), (n 1), ..., (n n)]
    */
    std::vector<double> generateBinomialCoefficients(int n);

}