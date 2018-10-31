include_directories(${PROJECT_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "-std=c++11")

add_executable(log log.cpp)
add_executable(single_layer_test single_layer_circle_test.cpp)

target_link_libraries(single_layer_test single_layer parametrizations quadrature)
target_link_libraries(log quadrature)
