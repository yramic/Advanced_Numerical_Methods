# Mandatory information that CMAKE needs
cmake_minimum_required(VERSION 2.6.0)

# Project name
project("cmTest")

# set cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/Modules/")

include("../cmake/Compiler.cmake")

#add eigen
include("../cmake/GetEigen.cmake")

#add Hilbert library
set( PROJECT_LINK_LIBS  libCppHilbert.so )
link_directories( ${CMAKE_BINARY_DIR}/../include )

set( CPPHILBERT_INCLUDE   ${CMAKE_SOURCE_DIR}/../Library )
include_directories(${CPPHILBERT_INCLUDE})

# Set up test for buildV, buildW, buildK, etc
set(SRC_test1 "test_build.cpp") # Source
add_executable(test1 ${SRC_test1})
add_dependencies(test1 Eigen)
target_link_libraries(test1 ${PROJECT_LINK_LIBS})

# Set up test for evaluateV, evaluateW, evaluateK, etc
set(SRC_test2 "test_evaluate.cpp") # Source
add_executable(test2 ${SRC_test2})
add_dependencies(test2 Eigen)
target_link_libraries(test2 ${PROJECT_LINK_LIBS})
