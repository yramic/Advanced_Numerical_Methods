cmake_minimum_required(VERSION 3.1)

project(LowRankApp)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(EIGEN3 REQUIRED)
if(NOT EIGEN3_FOUND)
	message(FATAL_ERROR "Could not find EIGEN3, consider setting EIGEN3_ROOT")
endif()

include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(include)
add_subdirectory(src)

add_executable(low_rank_test low_rank_test.cpp ${cpp_files} ${headers})

