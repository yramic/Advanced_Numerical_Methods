set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(${CMAKE_SOURCE_DIR}/include)

include_directories(../../CppHilbert/Library/source/)

include_directories(./googletest/include)

add_subdirectory(../../CppHilbert/Library ${CMAKE_BINARY_DIR}/CppHilbert)

link_directories(${CMAKE_BINARY_DIR}/CppHilbert)

# For using the g-test framework
add_subdirectory(googletest)

add_executable(tests tests.cpp)

target_link_libraries(tests parametrizations gtest single_layer quadrature CppHilbert)

#install(TARGETS param_tests_gtest DESTINATION ${CMAKE_SOURCE_DIR}/bin)

