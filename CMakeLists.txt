cmake_minimum_required(VERSION 3.0.0)
project(Bezier VERSION 0.1.0 DESCRIPTION "Library for Construction of Bezier Curves and Splines")

include_directories(${CMAKE_SOURCE_DIR}/include /usr/local/include)

add_library(Bezier SHARED include/bezier.h src/bezier_util.cpp ${CMAKE_SOURCE_DIR}/include/bezier_util.h)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (Bezier Eigen3::Eigen)

# Installing
set_target_properties(Bezier PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(Bezier PROPERTIES PUBLIC_HEADER include/bezier.h)
target_include_directories(Bezier PRIVATE .)
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS Bezier
        EXPORT BezierTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    BezierConfig.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
)
install(EXPORT BezierTargets
        FILE BezierConfig.cmake
        NAMESPACE Bezier::
        DESTINATION lib/cmake/Bezier)

add_library(nlopt SHARED IMPORTED)
set_target_properties(nlopt PROPERTIES IMPORTED_LOCATION /usr/local/lib/libnlopt.so)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(optimization src/optimization.cpp include/bezier.h src/bezier_util.cpp include/bezier_util.h)
target_link_libraries(optimization Eigen3::Eigen nlopt)

add_executable(visualization tests/visualization.cpp tests/matplotlibcpp.h include/bezier.h src/bezier_util.cpp include/bezier_util.h)
target_link_libraries(visualization Eigen3::Eigen ${PYTHON_LIBRARIES})

# Testing

include(CTest)
enable_testing()

add_executable(binomial_test tests/binomial_test.cpp include/bezier_util.h src/bezier_util.cpp)
add_executable(basic_curve_test include/bezier.h tests/basic_curve_test.cpp include/bezier_util.h src/bezier_util.cpp)
target_link_libraries (basic_curve_test Eigen3::Eigen)
add_executable(derivative_curve_test include/bezier.h tests/derivative_curve_test.cpp include/bezier_util.h src/bezier_util.cpp)
target_link_libraries (derivative_curve_test Eigen3::Eigen)
add_executable(interpolation_test tests/matplotlibcpp.h include/bezier.h tests/interpolation_test.cpp include/bezier_util.h src/bezier_util.cpp)
target_link_libraries (interpolation_test Eigen3::Eigen ${PYTHON_LIBRARIES})

add_test(binomial_test binomial_test COMMAND $<TARGET_FILE:cpp_test>)
add_test(basic_curve_test basic_curve_test COMMAND $<TARGET_FILE:cpp_test>)
add_test(derivative_curve_test derivative_curve_test COMMAND $<TARGET_FILE:cpp_test>)
add_test(interpolation_test interpolation_test COMMAND $<TARGET_FILE:cpp_test>)