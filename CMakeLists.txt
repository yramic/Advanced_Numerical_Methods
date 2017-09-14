cmake_minimum_required(VERSION 3.1)

project(AdvNumCSE)

###########################################################################
#
# INCLUDE MODULES AND FUNCTIONS FOR CMAKE
#
###########################################################################

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
include(${CMAKE_SOURCE_DIR}/cmake/functions/add_all_subdirectories.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/functions/add_executable_advnumcse.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/functions/get_target_name_advnumcse.cmake)

###########################################################################
#
# COMPILATION SETTINGS
#
###########################################################################

# DISALLOW IN-SOURCE BUILD
###########################################################################
#if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
#    message(FATAL_ERROR "Use an out of source build. Create a separate binary directory and run CMake there.")
#endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

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

include(${CMAKE_SOURCE_DIR}/cmake/functions/GetEIGEN3.cmake)

###########################################################################
#
# ADD SUBDIRECTORIES
#
###########################################################################

add_subdirectory(BEM)
add_subdirectory(HMAT)
add_subdirectory(CONVQUAD)

