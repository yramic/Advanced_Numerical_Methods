cmake_minimum_required(VERSION 3.10)

project(homeworks)  # only line that differs from developers/CMakeLists.txt

# Add custom targets
include("../cmake/modules/custom_targets.cmake")
get_custom_test_targets()

add_subdirectory(LowTriangToeplitz)
add_subdirectory(AbelIntegralEquation)
add_subdirectory(FractionalHeatEquation)
add_subdirectory(AbsorbingBoundaryCondition)
