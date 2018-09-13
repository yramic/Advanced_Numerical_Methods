include_directories(./../Library/src)

include_directories(./googletest/include)

set(CMAKE_CXX_FLAGS "-std=c++11")

# For using the g-test framework
add_subdirectory(googletest)

#add_executable(param_tests parametrization_tests.cpp)

#target_link_libraries(param_tests 2dparametricbem)

add_executable(param_tests_gtest parametrization_tests_gtest.cpp)

target_link_libraries(param_tests_gtest 2dparametricbem gtest)

