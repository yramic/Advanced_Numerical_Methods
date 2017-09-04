# Mandatory information that CMAKE needs
cmake_minimum_required( VERSION 3.0)

# Project name
project(Hilbert)

# Set Default build type
if( NOT CMAKE_BUILD_TYPE )
    message(STATUS "No build type selected. Default to DEBUG")
    set (CMAKE_BUILD_TYPE "DEBUG")
endif( )

# Set compiler flags
include("${CMAKE_SOURCE_DIR}/../cmake/Compiler.cmake")

# set cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/Modules/")

# Include Eigen
include("${CMAKE_SOURCE_DIR}/../cmake/GetEigen.cmake")

# Include CppHilbert
include_directories(source)

file(GLOB SOURCES "./*.hpp" "source/*.cpp" "source/*.hpp" "source/*.h" "source/*.c")

add_library(CppHilbert SHARED ${SOURCES})

# Define installation location
if (NOT CPPHILBERT_ROOT AND NOT $ENV{CPPHILBERT_ROOT} STREQUAL "")
  set( CPPHILBERT_LIB   $ENV{CPPHILBERT_ROOT}/include  )
else()
  set( CPPHILBERT_LIB   ${CMAKE_BINARY_DIR}/../include )
endif ()
install(TARGETS CppHilbert DESTINATION ${CPPHILBERT_LIB})



