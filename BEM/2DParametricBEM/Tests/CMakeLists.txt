include_directories(./../Library/src)

set(CMAKE_CXX_FLAGS "-std=c++11")

add_executable(param_tests parametrization_tests.cpp)

target_link_libraries(param_tests 2dparametricbem)

