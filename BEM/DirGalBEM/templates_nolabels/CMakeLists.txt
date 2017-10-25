cmake_minimum_required(VERSION 3.1)

project(DirichletBEM)

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
    add_custom_target(Eigen)
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
      SET(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)
endif()
add_definitions("-Wno-ignored-attributes -Wno-deprecated-declarations")
ExternalProject_Add(
    CppHilbertDownload
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../CppHilbert/Library
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/CppHilbert_install
    CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/CppHilbert_install -DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}/Eigen_install -DCMAKE_CXX_FLAGS=-I${EIGEN3_INCLUDE_DIR} # DONT DO THIS AT HOME
)
add_dependencies(CppHilbertDownload Eigen)
include_directories("../../CppHilbert/Library")
link_directories(${CMAKE_CURRENT_BINARY_DIR}/CppHilbertDownload-prefix/src/include/)
# make DBEMmain
set(SRCS DBEMmain.cpp)
add_executable(DBEMmain ${SRCS})
add_dependencies(DBEMmain Eigen CppHilbertDownload)
target_link_libraries(DBEMmain CppHilbert)
