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

set(SRC_test "test.cpp") # Source
add_executable(cmtest ${SRC_test})

add_dependencies(cmtest Eigen)
target_link_libraries(cmtest ${PROJECT_LINK_LIBS})
