set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(${CMAKE_SOURCE_DIR}/include)

FILE(GLOB ParamSources ./*.cpp)

add_library(parametrizations STATIC ${ParamSources})

install(TARGETS parametrizations DESTINATION ${CMAKE_SOURCE_DIR}/lib)

