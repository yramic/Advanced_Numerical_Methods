cmake_minimum_required(VERSION 3.10)

include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/cpp-pm/hunter/archive/v0.23.321.tar.gz"
    SHA1 "5e53cbb0429037ea8e2592bfd92704b8ff3ab492"
    LOCAL
)

project(ADVNUMCSE)

message("Project source dir =" "${PROJECT_SOURCE_DIR}")
message("Project binary dir =" "${PROJECT_BINARY_DIR}")

# Check that this is not an in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

#setup compiler
set(CMAKE_CXX_STANDARD 17)

#enable warnings
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
CHECK_CXX_COMPILER_FLAG(-Wextra COMPILER_SUPPORTS_WEXTRA)
CHECK_CXX_COMPILER_FLAG(-Wconversion COMPILER_SUPPORTS_WCONVERSION)

add_library(warning_flags INTERFACE )
if(COMPILER_SUPPORTS_WALL)
  target_compile_options(warning_flags INTERFACE -Wall)
endif()
if(COMPILER_SUPPORTS_WEXTRA)
  target_compile_options(warning_flags INTERFACE -Wextra)
endif()
if(COMPILER_SUPPORTS_WCONVERSION)
  target_compile_options(warning_flags INTERFACE -Wconversion)
endif()

# Add switches
###############################################################################
option(LEHRFEM "usage: cmake -DLEHRFEM=ON/OFF .." OFF)
message(STATUS "LEHRFEM ${LEHRFEM}")

option(LECTURECODES "usage: cmake -DLECTURECODES=ON/OFF .." ON)
message(STATUS "LECTURECODES ${LECTURECODES}")

option(MASTERSOLUTION "usage: cmake -DMASTERSOLUTION=ON/OFF .." ON)
message(STATUS "MASTERSOLUTION ${MASTERSOLUTION}")

option(MYSOLUTION "usage: cmake -DMYSOLUTION=ON/OFF .." ON)
message(STATUS "MYSOLUTION ${MYSOLUTION}")

# Get Dependencies
# (don't forget to update cmake/Config.cmake.in !
###############################################################################
# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)

# Get LehrFEMpp
if(LEHRFEM)
  hunter_add_package(lehrfempp)
  find_package(lehrfempp CONFIG REQUIRED)
endif()
# Get Google Test
hunter_add_package(GTest)
find_package(GTest CONFIG REQUIRED)

# Get some boost program options if LectureCodes are enabled:
if(LECTURECODES)
  hunter_add_package(Boost COMPONENTS program_options)
  find_package(Boost CONFIG REQUIRED program_options)
endif()

# Add subdirectories
if(LEHRFEM)
  add_subdirectory(DUMMYCHAPTER)
endif()
add_subdirectory(HMAT)
add_subdirectory(CONVQUAD)
add_subdirectory(AMG)
