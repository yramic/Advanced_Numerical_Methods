cmake_minimum_required(VERSION 3.1)

project(HMAT_1D)

add_subdirectory(include)
add_subdirectory(src)

add_executable(low_rank_test_1D
               low_rank_test.cpp
               ${cpp_files} ${headers})

add_executable(low_rank_test_1D_uni
               low_rank_test_uni.cpp
               ${cpp_files_uni} ${cpp_files}
               ${headers_uni} ${headers})
