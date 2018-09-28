set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(./BEM\ Space)

include_directories(Quadrature)

include_directories(../../../CppHilbert/Library/source/)

#External library for Gauss Quadrature
#include_directories(/home/piyush/Downloads/cpp/src)
#FILE(GLOB ExtSources /home/piyush/Downloads/cpp/src/*.cpp)

add_library(parametrizations STATIC parametrized_fourier_sum.cpp parametrized_line.cpp parametrized_semi_circle.cpp)

add_library(single_layer STATIC single_layer.cpp ./Quadrature/logweight_quadrature.cpp ./Quadrature/genLaguerreRule.cpp ../../../CppHilbert/Library/source/gaussQuadrature.c)
