#include "util.h"

namespace Util {

    std::map<int, std::vector<double>> pascalTriangleMap; 

    /**
     * https://stackoverflow.com/questions/15580291/how-to-efficiently-calculate-a-row-in-pascals-triangle
    */
    std::vector<double> generateBinomialRow(int n) {
        std::vector<double> row = {1};
        for (int k=0; k < n; k++) {
            row.push_back(row[k] * (n-k) / (k+1));
        }
        return row;
    }   

    std::vector<double> generateBinomialCoefficients(int n) {
        if (!pascalTriangleMap.count(n)) {
            pascalTriangleMap[n] = generateBinomialRow(n);
        }
        return pascalTriangleMap[n];
    }

}