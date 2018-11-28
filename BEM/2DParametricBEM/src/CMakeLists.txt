set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(Parametrizations)

add_subdirectory(Quadrature)

add_library(single_layer STATIC single_layer.cpp parametrized_mesh.cpp)
add_library(double_layer STATIC double_layer.cpp parametrized_mesh.cpp)
add_library(hypersingular STATIC hypersingular.cpp parametrized_mesh.cpp)
