cmake_minimum_required(VERSION 3.1)

project(Hilbert)

include_directories(source)

file(GLOB SOURCES "./*.hpp" "source/*.cpp" "source/*.hpp" "source/gaussTable.h" "source/gaussQuadrature.h" "source/gaussQuadrature.c")

add_library(CppHilbert SHARED ${SOURCES})

# Define installation location
if(NOT CPPHILBERT_ROOT AND NOT $ENV{CPPHILBERT_ROOT} STREQUAL "")
    set( CPPHILBERT_LIB   $ENV{CPPHILBERT_ROOT}/include  )
else()
    set( CPPHILBERT_LIB   ${CMAKE_BINARY_DIR}/../include )
endif()
install(TARGETS CppHilbert DESTINATION ${CPPHILBERT_LIB})

