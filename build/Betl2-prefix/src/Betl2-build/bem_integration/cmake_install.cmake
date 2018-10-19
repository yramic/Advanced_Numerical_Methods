# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libbem_integration.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_integration/libbem_integration.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/cache.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/enum_singularity.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/fix_edge_orientations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_base_integrator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_integrator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_integrator_functors.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_integrator_regular.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_integrator_singular.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel_evaluator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel_layer_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel_layer_traits_edge.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel_layer_traits_lagrange.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel_layer_traits_mixed.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_kernel_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_local_evaluations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_local_evaluations_lagrange.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_local_evaluations_mixed.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_quadrature.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_quadrature_creators.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_quadrature_rules.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_singular_storage.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_singularity_detector.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_static_cache.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/galerkin_static_cache_evaluator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/integration_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/sautertransformation.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/sautertransformation_quad_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/sautertransformation_quad_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/sautertransformation_tria_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/sautertransformation_tria_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/singular_quadrature_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration/singularity_mapping.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/bem_integration" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/cache.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/enum_singularity.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/fix_edge_orientations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_base_integrator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_integrator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_integrator_functors.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_integrator_regular.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_integrator_singular.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel_evaluator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel_layer_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel_layer_traits_edge.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel_layer_traits_lagrange.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel_layer_traits_mixed.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_kernel_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_local_evaluations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_local_evaluations_lagrange.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_local_evaluations_mixed.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_quadrature.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_quadrature_creators.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_quadrature_rules.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_singular_storage.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_singularity_detector.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_static_cache.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/galerkin_static_cache_evaluator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/integration_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/sautertransformation.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/sautertransformation_quad_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/sautertransformation_quad_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/sautertransformation_tria_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/sautertransformation_tria_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/singular_quadrature_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_integration/singularity_mapping.hpp"
    )
endif()

