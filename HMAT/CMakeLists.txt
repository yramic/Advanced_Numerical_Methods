cmake_minimum_required(VERSION 3.10)

project(homeworks)  # only line that differs from developers/CMakeLists.txt

# Add custom targets
include("../cmake/modules/custom_targets.cmake")
get_custom_test_targets()

add_subdirectory(CLUSTERING) # Builds a library
add_subdirectory(CLUSTERING_revised) # Builds a library
add_subdirectory(KernMatLLRApprox) # Uses a library build in CLUSTERING
add_subdirectory(KernMatLLRApprox_revised) # Uses a library build in CLUSTERING

add_subdirectory(TensorProductChebIntp)
add_subdirectory(LowRankMerge)
add_subdirectory(GravitationalForces)
add_subdirectory(EvalTrigPoly)
