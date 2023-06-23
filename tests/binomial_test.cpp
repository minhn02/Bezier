#include <cassert>
#include <vector>
#include <cstdio>
#include <chrono>

#include "bezier.h"

std::vector<std::vector<double>> generatePascalTriangleToRow10() {
    std::vector<std::vector<double>> triangle;
    triangle.push_back({1});
    triangle.push_back({1, 1});
    triangle.push_back({1, 2, 1});
    triangle.push_back({1, 3, 3, 1});
    triangle.push_back({1, 4, 6, 4, 1});
    triangle.push_back({1, 5, 10, 10, 5, 1});
    triangle.push_back({1, 6, 15, 20, 15, 6, 1});
    triangle.push_back({1, 7, 21, 35, 35, 21, 7, 1});
    triangle.push_back({1, 8, 28, 56, 70, 56, 28, 8, 1});
    triangle.push_back({1, 9, 36, 84, 126, 126, 84, 36, 9, 1});
    triangle.push_back({1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1});
    
    return triangle;
}

void test_binomial_generation() {
    std::vector<std::vector<double>> triangle = generatePascalTriangleToRow10();
    for (int i=0; i < 10; i++) {
        assert(triangle[i] == Bezier::generateBinomialCoefficients(i));
    }
}

int main(int argc, char const *argv[]) {
    auto before = std::chrono::high_resolution_clock::now();
    test_binomial_generation();
    auto after = std::chrono::high_resolution_clock::now();
    auto non_cache_time = std::chrono::duration_cast<std::chrono::nanoseconds>(after-before);

    // check caching
    before = std::chrono::high_resolution_clock::now();
    test_binomial_generation();
    after = std::chrono::high_resolution_clock::now();
    auto cache_time = std::chrono::duration_cast<std::chrono::nanoseconds>(after-before);

    // cached time should be less than initial call
    assert(cache_time.count() < non_cache_time.count());
}
