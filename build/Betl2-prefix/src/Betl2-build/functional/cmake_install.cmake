# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libfunctional.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/functional/libfunctional.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/L_two_product_evaluator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/analytical_grid_function.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/bem_analytical_grid_function.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/dof_interpolator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/dof_on_element_interpolator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/dof_on_element_interpolator_curl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/dof_on_element_interpolator_div.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/dof_on_element_interpolator_lagrange.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/dof_on_element_interpolator_lagrange_hierarchical.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/edge_gridfunction_wrapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/edge_integrator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/edge_moments.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/edge_parametrizations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/face_integrator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/face_moments.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/face_parametrizations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/fem_analytical_grid_function.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/grid_function.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/grid_function_integrator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/grid_function_operations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/grid_function_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/h1_moments.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/interpolate_edge_dofs.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/interpolate_h1_dofs.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/interpolation_grid_function.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/interpolation_grid_function_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/interpolation_kernels.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/representation_solution.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/trace_enumerators.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/unisolvent_interpolation_curl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/unisolvent_interpolation_div.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/unisolvent_interpolation_lagrange_hierarchical.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/unisolvent_interpolation_scheme.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/volume_integrator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/functional/volume_moments.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/functional" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/L_two_product_evaluator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/analytical_grid_function.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/bem_analytical_grid_function.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/dof_interpolator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/dof_on_element_interpolator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/dof_on_element_interpolator_curl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/dof_on_element_interpolator_div.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/dof_on_element_interpolator_lagrange.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/dof_on_element_interpolator_lagrange_hierarchical.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/edge_gridfunction_wrapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/edge_integrator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/edge_moments.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/edge_parametrizations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/face_integrator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/face_moments.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/face_parametrizations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/fem_analytical_grid_function.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/grid_function.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/grid_function_integrator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/grid_function_operations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/grid_function_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/h1_moments.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/interpolate_edge_dofs.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/interpolate_h1_dofs.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/interpolation_grid_function.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/interpolation_grid_function_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/interpolation_kernels.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/representation_solution.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/trace_enumerators.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/unisolvent_interpolation_curl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/unisolvent_interpolation_div.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/unisolvent_interpolation_lagrange_hierarchical.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/unisolvent_interpolation_scheme.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/volume_integrator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/functional/volume_moments.hpp"
    )
endif()

