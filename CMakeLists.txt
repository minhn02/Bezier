cmake_minimum_required(VERSION 3.0.0)
project(Bezier VERSION 0.1.0)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(Bezier src/bezier.cpp include/bezier.h src/util.cpp src/util.h)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (Bezier Eigen3::Eigen)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Testing

include(CTest)
enable_testing()

add_executable(binomial_test tests/binomial_test.cpp src/util.h src/util.cpp)
add_executable(basic_curve_test include/bezier.h src/bezier.cpp tests/basic_curve_test.cpp src/util.h src/util.cpp)
target_link_libraries (basic_curve_test Eigen3::Eigen)
add_executable(derivative_curve_test include/bezier.h src/bezier.cpp tests/derivative_curve_test.cpp src/util.h src/util.cpp)
target_link_libraries (derivative_curve_test Eigen3::Eigen)

add_test(binomial_test binomial_test COMMAND $<TARGET_FILE:cpp_test>)
add_test(basic_curve_test basic_curve_test COMMAND $<TARGET_FILE:cpp_test>)
add_test(derivative_curve_test derivative_curve_test COMMAND $<TARGET_FILE:cpp_test>)