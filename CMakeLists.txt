cmake_minimum_required(VERSION 3.1)

project(AdvNumCSE)

###########################################################################
#
# INCLUDE MODULES AND FUNCTIONS FOR CMAKE
#
###########################################################################

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/modules")
include(${CMAKE_SOURCE_DIR}/CMake/functions/add_all_subdirectories.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/functions/add_executable_advnumcse.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/functions/add_executable_advnumcse_betl2.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/functions/get_target_name_advnumcse.cmake)

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

# STORE BINARIES IN ./bin
#(else CMake cannot find files to link for automatic targets)
###########################################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# SET DEFAULT BUILD TYPE
###########################################################################
if( NOT CMAKE_BUILD_TYPE )
    message(STATUS "No build type selected. Default to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# COMPILATION FLAGS
###########################################################################
set(CMAKE_CXX_STANDARD 11)

# IGNORE WARNINGS
###########################################################################

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_definitions(-Wno-deprecated-register -Wno-return-type-c-linkage -Wno-keyword-macro -Wno-missing-braces)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    add_definitions(-Wno-deprecated-declarations -Wno-ignored-attributes -Wno-misleading-indentation -Wno-unknown-pragmas)
endif()

###########################################################################
#
# FIND DEPENDENCIES
#
###########################################################################

# FIND EIGEN
###########################################################################

include(${CMAKE_SOURCE_DIR}/CMake/functions/GetEIGEN3.cmake)
add_definitions(-DEIGEN3_ACTIVATED)

###########################################################################
#
# INCLUDE DIRECTORIES
#
###########################################################################
include_directories("BEM/CppHilbert/Library")
include_directories("third_party")

###########################################################################
#
# BETL AND ETHGENERICGRID
#
###########################################################################

# SAFETY CHECKUPS
###########################################################################

# There is too much hassle handling these corner cases, so we just require the
# people who have Betl installed to undefine these
if (DEFINED ENV{ETH_ROOT})
  message(FATAL_ERROR "EE: Environment variable ETH_ROOT found. Please unset using\n    unset ETH_ROOT")
endif()

if (DEFINED ENV{BETL2_ROOT})
  message(FATAL_ERROR "EE: Environment variable BETL2_ROOT found. Please unset using\n    unset BETL2_ROOT")
endif()

# Make sure they have checked out the repository
if (NOT EXISTS ${CMAKE_SOURCE_DIR}/third_party/Betl2/CMakeLists.txt)
  message(FATAL_ERROR "You have not downloaded the Betl2 submodule. Please do\n    cd ${CMAKE_SOURCE_DIR}/third_party/\n    git checkout Betl2;\n    cd Betl2;\n    git submodule init;\n    git submodule update --init --recursive\n\nRemember you have to have a gitlab.math.ethz.ch account for this. Go to https://www.gitlab.math.ethz.ch and log in to get an account, then send an email to the organizers to get access.")
endif()

if (NOT EXISTS ${CMAKE_SOURCE_DIR}/third_party/Betl2/3rd_Party/ethGenericGrid/CMakeLists.txt)
  message(FATAL_ERROR "You have not downloaded the ethGenericGrid submodule. Please do\n    cd ${CMAKE_SOURCE_DIR}/third_party/Betl2/3rd_Party\n    git checkout ethGenericGrid;\n    cd ethGenericGrid;\n    git submodule init;\n    git submodule update --init --recursive\nRemember you have to have a gitlab.math.ethz.ch account for this. Go to https://www.gitlab.math.ethz.ch and log in to get an account, then send an email to the organizers to get access.")
endif()

set(ETH_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install/include)
set(ETH_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install/lib)

# ETH GENERIC GRID
###########################################################################
ExternalProject_Add(
    ethGenericGrid
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third_party/Betl2/3rd_Party/ethGenericGrid
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install
    CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install -DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}/Eigen_install -DETH_ROOT=${CMAKE_CURRENT_BINARY_DIR}/ethGenericGrid_install -DCMAKE_INCLUDE_DIRECTORIES_BEFORE=${NUMPDE_INCLUDE_DIRECTORIES_BEFORE}
)

add_dependencies(ethGenericGrid Eigen)

# BETL2
###########################################################################
ExternalProject_Add(
    Betl2
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third_party/Betl2/Library
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

find_package(Boost REQUIRED COMPONENTS program_options)
if (${Boost_MAJOR_VERSION} VERSION_LESS 1 OR ${Boost_MINOR_VERSION} VERSION_LESS 58)
  message(STATUS "Downloading odeint, since boost version too old.")
  ExternalProject_Add(
    ODEINT
    URL https://github.com/headmyshoulder/odeint-v2/archive/v2.6_boost1.60.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/odeint_source
    INSTALL_COMMAND  ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    )
  include_directories (${CMAKE_CURRENT_BINARY_DIR}/odeint_source/include)
else()
  add_custom_target(ODEINT) # dependency dummy
endif()
include_directories (${Boost_INCLUDE_DIRS})

set(BETL2_LIBRARIES bem_integration
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
${Boost_LIBRARIES}
)

###########################################################################
#
# ADD SUBDIRECTORIES
#
###########################################################################

add_subdirectory(BEM)
add_subdirectory(HMAT)
add_subdirectory(CONVQUAD)

