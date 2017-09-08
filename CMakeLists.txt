cmake_minimum_required(VERSION 3.1)

project(AdvNumCSE)

###########################################################################
#
# INCLUDE MODULES AND FUNCTIONS FOR CMAKE
#
###########################################################################

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

###########################################################################
#
# COMPILATION SETTINGS
#
###########################################################################

# DISALLOW IN-SOURCE BUILD
###########################################################################
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "Use an out of source build. Create a separate binary directory and run CMake there.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# SET DEFAULT BUILD TYPE
###########################################################################
if( NOT CMAKE_BUILD_TYPE )
    message(STATUS "No build type selected. Default to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# COMPILATION FLAGS
###########################################################################
set(CMAKE_CXX_STANDARD 11)

###########################################################################
#
# FIND DEPENDENCIES
#
###########################################################################

# FIND EIGEN
###########################################################################

include(ExternalProject) # To download Eigen if it is missing

find_package(EIGEN3)

if(${EIGEN3_FOUND})

    include_directories(${EIGEN3_INCLUDE_DIR})

else()

    SET(DOWNLOADING_EIGEN ON)
    message("-- Downloading Eigen3")
    ExternalProject_Add(
        Eigen
        URL http://bitbucket.org/eigen/eigen/get/3.3.4.zip
         SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
        DOWNLOAD_NO_PROGRESS 1
        CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install)
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)

endif()

###########################################################################
#
# ADD SUBDIRECTORIES
#
###########################################################################

add_subdirectory(LowRank)

