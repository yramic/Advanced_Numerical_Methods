# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/excalibur/AdvNum/Code/build/betl2_install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DEBUG")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/direct_solver.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/identity_matrix.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/lapack_declarations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/lapack_defines.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/lapack_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/linalg_settings.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/no_preconditioner.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/permutation_wrapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparse_assign_policy.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparse_assign_policy_petsc.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparse_enumerators.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparsity_manager.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparsity_options.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparsity_pattern.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparsity_policy.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/sparsity_policy_eigen.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/triplet_operators.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg/very_sparse_matrix.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/linalg" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/direct_solver.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/identity_matrix.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/lapack_declarations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/lapack_defines.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/lapack_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/linalg_settings.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/no_preconditioner.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/permutation_wrapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparse_assign_policy.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparse_assign_policy_petsc.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparse_enumerators.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparsity_manager.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparsity_options.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparsity_pattern.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparsity_policy.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/sparsity_policy_eigen.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/triplet_operators.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/linalg/very_sparse_matrix.hpp"
    )
endif()

