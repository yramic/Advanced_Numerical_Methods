include_directories(${CMAKE_SOURCE_DIR}/include)

include_directories(../../../CppHilbert/Library/source/)

FILE(GLOB QuadratureSources ./*.cpp)

add_library(quadrature STATIC ${QuadratureSources})



