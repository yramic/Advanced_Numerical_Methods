cmake_minimum_required(VERSION 3.1)

project(cmTest)

set(SRC_test
    test_build.cpp)


add_executable(cmtest ${SRC_test})
#add_dependencies(cmtest Eigen) # already included in top-level CMakeLists.txt
target_link_libraries(cmtest CppHilbert)

