cmake_minimum_required(VERSION 3.1)

project(HMAT_1D)

add_subdirectory(include)
add_subdirectory(src)

add_library(hmat_1d SHARED ${headers} ${cpp_files})
add_executable(low_rank_test_1D low_rank_test.cpp)
target_link_libraries(low_rank_test_1D hmat_1d)

add_library(hmat_1d_uni_sols SHARED ${headers} ${headers_uni} ${cpp_files_uni} ${cpp_files_uni_sols})
add_executable(low_rank_test_1D_uni_solution low_rank_test_uni.cpp)
target_link_libraries(low_rank_test_1D_uni_solution hmat_1d hmat_1d_uni_sols)

