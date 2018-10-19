# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/combinatorial_divergence.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/combinatorial_divergence_element_matrices.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/combinatorial_divergence_sp.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/combinatorial_gradient.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/combinatorial_gradient_element_matrices.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/combinatorial_gradient_sp.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/identity_operator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/identity_operator_sp.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/identity_operator_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/p_restriction_operator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/p_restriction_operator_sp.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/sparse_base_operator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/surface_restriction_edge_mapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/surface_restriction_operator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/surface_restriction_operator_sp.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/vector_basis_grid_function.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators/vector_h1_interpolator_sp.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/sparse_operators" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/combinatorial_divergence.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/combinatorial_divergence_element_matrices.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/combinatorial_divergence_sp.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/combinatorial_gradient.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/combinatorial_gradient_element_matrices.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/combinatorial_gradient_sp.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/identity_operator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/identity_operator_sp.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/identity_operator_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/p_restriction_operator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/p_restriction_operator_sp.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/sparse_base_operator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/surface_restriction_edge_mapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/surface_restriction_operator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/surface_restriction_operator_sp.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/vector_basis_grid_function.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/sparse_operators/vector_h1_interpolator_sp.hpp"
    )
endif()

