set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(./../Library/src)

include_directories(./../Library/src/BEM\ Space)

include_directories(./googletest/include)

add_subdirectory(../../CppHilbert/Library ${CMAKE_SOURCE_DIR}/build/bin)

link_directories(${CMAKE_SOURCE_DIR}/build/bin)

include_directories(../../CppHilbert/Library/source/)

# For using the g-test framework
add_subdirectory(googletest)

add_executable(param_tests_gtest parametrization_tests_gtest.cpp)

target_link_libraries(param_tests_gtest parametrizations gtest single_layer CppHilbert)

