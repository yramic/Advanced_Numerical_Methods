cmake_minimum_required( VERSION 2.8 )

set(CMAKE_CXX_FLAGS "-std=c++11")

add_library(2dparametricbem STATIC parametrized_fourier_sum.cpp parametrized_line.cpp parametrized_semi_circle.cpp)

