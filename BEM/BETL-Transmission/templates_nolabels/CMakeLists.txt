cmake_minimum_required(VERSION 3.1)

project(BETL-Debug)

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

      endif()

      ADD_DEFINITIONS(-DEIGEN3_ACTIVATED)

add_definitions("-Wno-ignored-attributes -Wno-deprecated-declarations")


# There is too much hassle handling these corner cases, so we just require the people who have
# Betl installed to undefine these
if (DEFINED ENV{ETH_ROOT})
  message(FATAL_ERROR "EE: Environment variable ETH_ROOT found. Please unset using\n    unset ETH_ROOT")
endif()

if (DEFINED ENV{BETL2_ROOT})
  message(FATAL_ERROR "EE: Environment variable BETL2_ROOT found. Please unset using\n    unset BETL2_ROOT")
endif()

# Make sure they have checked out the repository

if (NOT EXISTS ${CMAKE_SOURCE_DIR}/../../../third_party/Betl2/CMakeLists.txt)
  message(FATAL_ERROR "You have not downloaded the Betl2 submodule. Please do\n    cd ${CMAKE_SOURCE_DIR}/third_party/\n    git checkout Betl2;\n    cd Betl2;\n    git submodule init;\n    git submodule update --init --recursive\n\nRemember you have to have a gitlab.math.ethz.ch account for this. Go to https://www.gitlab.math.ethz.ch and log in to get an account, then send an email to the organizers to get access.")
endif()


if (NOT EXISTS ${CMAKE_SOURCE_DIR}/../../../third_party/Betl2/3rd_Party/ethGenericGrid/CMakeLists.txt)
  message(FATAL_ERROR "You have not downloaded the ethGenericGrid submodule. Please do\n    cd ${CMAKE_SOURCE_DIR}/third_party/Betl2/3rd_Party\n    git checkout ethGenericGrid;\n    cd ethGenericGrid;\n    git submodule init;\n    git submodule update --init --recursive\nRemember you have to have a gitlab.math.ethz.ch account for this. Go to https://www.gitlab.math.ethz.ch and log in to get an account, then send an email to the organizers to get access.")
endif()

set(ETH_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install/include)
set(ETH_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install/lib)

ExternalProject_Add(
    ethGenericGrid
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../../third_party/Betl2/3rd_Party/ethGenericGrid
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install
    CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install -DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}/Eigen_install -DETH_ROOT=${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install -DCMAKE_INCLUDE_DIRECTORIES_BEFORE=${NUMPDE_INCLUDE_DIRECTORIES_BEFORE}
)


add_dependencies(ethGenericGrid Eigen)


# include betl
ExternalProject_Add(
    Betl2
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../../third_party/Betl2/Library
    CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/betl2_install
    -DETH_INCLUDE_DIRS=${ETH_INCLUDE_DIRS}/eth ${CMAKE_BINARY_DIR}/Eigen_install/include -DETH_BASE_LIB=${ETH_LIBRARY_DIR} -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/betl2_install -DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}/Eigen_install -DETH_ROOT=${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install
-DCMAKE_INCLUDE_DIRECTORIES_BEFORE=${NUMPDE_INCLUDE_DIRECTORIES_BEFORE}
    )
add_dependencies(Betl2 ethGenericGrid)
add_dependencies(Betl2 Eigen)

include_directories(${ETH_INCLUDE_DIRS} ${ETH_INCLUDE_DIRS}/eth)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/betl2_install/include)
link_directories(${ETH_LIBRARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR}/betl2_install/lib)


SET(BETL2_LIBRARIES bem_integration
cmdl_parser
fe
functional
fundsol
gmsh_input
grid
hypre_bindings
material
petsc_bindings
utils
vtu_exporter
eth_eth_base
eth_interface
${ETH_LIBRARIES} 
)

set(SRCS main.cpp)
add_executable(main ${SRCS})
target_link_libraries(main ${BETL2_LIBRARIES})
add_dependencies(main Eigen)
add_dependencies(main Betl2)

