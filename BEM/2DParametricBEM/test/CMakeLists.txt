set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(../../CppHilbert/Library/source/)
include_directories(./googletest/include)

add_subdirectory(../../CppHilbert/Library ${CMAKE_SOURCE_DIR}/lib)

link_directories(${CMAKE_SOURCE_DIR}/lib)

# For using the g-test framework
add_subdirectory(googletest)

add_executable(param_tests_gtest parametrization_tests_gtest.cpp)

target_link_libraries(param_tests_gtest parametrizations gtest single_layer CppHilbert)

