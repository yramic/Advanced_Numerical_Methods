cmake_minimum_required(VERSION 3.1)

project(LogWeightedGaussRule)

###########################################################################
#
# INCLUDE MODULES AND FUNCTIONS FOR CMAKE
#
###########################################################################

# Notice double brackets for '${CMAKE_SOURCE_DIR}':
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

###########################################################################
#
# COMPILATION SETTINGS
#
###########################################################################

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

# Notice double brackets for '${EIGEN3_FOUND}':
if(${EIGEN3_FOUND})

    # Notice double brackets for '${EIGEN3_INCLUDE_DIR}':
    include_directories(${EIGEN3_INCLUDE_DIR})

else()

    # Notice double brackets for '${}':
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


# make main
FILE(GLOB SRC *.cpp *.hpp)
add_executable(main ${SRC})

