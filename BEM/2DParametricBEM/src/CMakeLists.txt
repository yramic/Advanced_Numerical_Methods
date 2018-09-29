set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(${CMAKE_SOURCE_DIR}/include)

include_directories(../../CppHilbert/Library/source/)

add_subdirectory(Parametrizations)

#External library for Gauss Quadrature
#include_directories(/home/piyush/Downloads/cpp/src)
#FILE(GLOB ExtSources /home/piyush/Downloads/cpp/src/*.cpp)

FILE(GLOB QuadratureSources ./Quadrature/*.cpp)

add_library(single_layer STATIC single_layer.cpp ${QuadratureSources} ../../CppHilbert/Library/source/gaussQuadrature.c)

#install(TARGETS single_layer DESTINATION ${CMAKE_SOURCE_DIR}/lib)
