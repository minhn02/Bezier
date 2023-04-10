cmake_minimum_required(VERSION 3.0.0)
project(Bezier VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(Bezier src/bezier.cpp include/bezier.h)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 
target_link_libraries (Bezier Eigen3::Eigen)