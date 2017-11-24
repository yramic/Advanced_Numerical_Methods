cmake_minimum_required(VERSION 3.1)

project(HMAT_BEM)

add_subdirectory(include)
add_subdirectory(src)

# Main target to get familiar with HMAT_BEM code
# Template target
add_library(hmat_bem_tmps SHARED ${headers} ${cpp_files} ${cpp_files_tmps})
target_link_libraries(hmat_bem_tmps CppHilbert)
add_executable(low_rank_test_BEM_template low_rank_test.cpp)
target_link_libraries(low_rank_test_BEM_template hmat_bem_tmps)

# Solution target
add_library(hmat_bem_sols SHARED ${headers} ${cpp_files} ${cpp_files_sols})
target_link_libraries(hmat_bem_sols CppHilbert)
add_executable(low_rank_test_BEM_solution low_rank_test.cpp)
target_link_libraries(low_rank_test_BEM_solution hmat_bem_sols)

