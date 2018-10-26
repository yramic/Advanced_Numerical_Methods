include_directories(${CMAKE_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "-std=c++11")

add_executable(log log.cpp)

target_link_libraries(log quadrature)
