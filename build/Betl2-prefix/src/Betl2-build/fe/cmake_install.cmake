# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libfe.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/fe/libfe.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_constraint.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_handler.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_data_set_factories.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_distribution_policy.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_distribution_policy_discontinuous.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_distribution_policy_continuous.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/dof_iterator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/domain_dof_marker.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/elemwise_constrained_fespace.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/filter_dofs.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/boundary_dof_marker.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/intersections_dof_marker.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/index_pair.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_enumerators.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/febasis.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fespace.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fespace_orientation_builder.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/condition_dof.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/constrained_fespace.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/febasis_utilities.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/febasis_to_reference_points.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/febasis_to_reference_points_util.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_basis_type_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_segment.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_basis_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_segment.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_lagrange_hierarchical_basis_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_basis_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_hierarchical_basis_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_hierarchical_basis_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_hierarchical_basis_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_hierarchical_basis_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_hierarchical_basis_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_gradient_lagrange_hierarchical_basis_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/MACRO_make_lagrange_basis.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_basis_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_curl_basis_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_curl_basis_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_curl_basis_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_curl_basis_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_curl_basis_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_curl_basis_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_curl_lagrange_basis_functions.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_linear_combination.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_linear_combination_lagrange.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_linear_combination_lagrange_hierarchical.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_linear_combination_curl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_linear_combination_div.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_linear_combination_higher_order_helper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/fe_multiplicity.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/hierarchical_febasis_mapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/impose_linear_combinations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/surface_dof_marker.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/fe/special_surface_dof_marker.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/fe" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_constraint.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_handler.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_data_set_factories.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_distribution_policy.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_distribution_policy_discontinuous.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_distribution_policy_continuous.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/dof_iterator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/domain_dof_marker.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/elemwise_constrained_fespace.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/filter_dofs.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/boundary_dof_marker.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/intersections_dof_marker.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/index_pair.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_enumerators.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/febasis.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fespace.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fespace_orientation_builder.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/condition_dof.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/constrained_fespace.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/febasis_utilities.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/febasis_to_reference_points.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/febasis_to_reference_points_util.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_basis_type_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_segment.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_basis_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_segment.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_lagrange_hierarchical_basis_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_basis_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_hierarchical_basis_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_hierarchical_basis_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_hierarchical_basis_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_hierarchical_basis_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_hierarchical_basis_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_gradient_lagrange_hierarchical_basis_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/MACRO_make_lagrange_basis.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_basis_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_curl_basis_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_curl_basis_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_curl_basis_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_curl_basis_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_curl_basis_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_curl_basis_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_curl_lagrange_basis_functions.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_linear_combination.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_linear_combination_lagrange.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_linear_combination_lagrange_hierarchical.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_linear_combination_curl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_linear_combination_div.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_linear_combination_higher_order_helper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/fe_multiplicity.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/hierarchical_febasis_mapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/impose_linear_combinations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/surface_dof_marker.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/fe/special_surface_dof_marker.hpp"
    )
endif()

