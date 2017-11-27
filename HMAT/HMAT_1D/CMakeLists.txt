cmake_minimum_required(VERSION 3.1)

project(HMAT_1D)

add_subdirectory(include)
add_subdirectory(src)

# Main target to get familiar with HMAT_1D code
# Template target for exercises e and g
add_library(hmat_1d_tmps SHARED ${headers} ${cpp_files} ${cpp_files_tmps})
add_executable(low_rank_test_1D_template low_rank_test.cpp)
target_link_libraries(low_rank_test_1D_template hmat_1d_tmps)

# Solution target for exercises e and g
add_library(hmat_1d_sols SHARED ${headers} ${cpp_files} ${cpp_files_sols})
add_executable(low_rank_test_1D_solution low_rank_test.cpp)
target_link_libraries(low_rank_test_1D_solution hmat_1d_sols)

# Template target for exercise b
add_library(hmat_1d_uni_tmps SHARED ${headers} ${headers_uni}
    ${cpp_files} ${cpp_files_tmps} ${cpp_files_uni} ${cpp_files_uni_tmps})
add_executable(low_rank_test_1D_uni_template low_rank_test_uni.cpp)
target_link_libraries(low_rank_test_1D_uni_template hmat_1d_uni_tmps)

# Solution target for exercise b
add_library(hmat_1d_uni_sols SHARED ${headers} ${headers_uni}
    ${cpp_files} ${cpp_files_sols} ${cpp_files_uni} ${cpp_files_uni_sols})
add_executable(low_rank_test_1D_uni_solution low_rank_test_uni.cpp)
target_link_libraries(low_rank_test_1D_uni_solution hmat_1d_uni_sols)

