# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libgrid.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/grid/libgrid.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/MACRO_grid_view.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/MACRO_grid_view_refined.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/bary_entities.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/bary_entity_element.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/bary_grid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/bary_grid_view.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/bary_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/bisection.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/edge_creator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/edges_on_face.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/element_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/element_typeinfo.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/element_typelist.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/element_types.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entities.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entity_edge.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entity_element.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entity_face.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entity_iterator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entity_node.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/entity_pointer.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/generic_coordinate_transformator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/geometry_identity_mapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/geometry_mapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/geometry_sphere_mapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entities.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entity_edge.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entity_element.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entity_iterator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entity_node.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entity_pointer.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid2d_entity_seed.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid_creator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid_partitioner.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/grid_view.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/index_set.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_iterator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_base_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_mapping_data.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/intersection_set.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/make_entity_pair.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/mapper_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/no_intersection_filter.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/num_nodes_to_grid2d_element_type.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/orientation_container.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/orientation_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/permutation.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/surface_grid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/surface_grid_init_from_input.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/surface_grid_init_from_other.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/surface_grid_init_from_volume_grid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/surface_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/volume2d_grid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/volume2d_grid_init_from_input.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/volume2d_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/volume_to_surface_mapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/grid/volume_traits.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/grid" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/MACRO_grid_view.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/MACRO_grid_view_refined.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/bary_entities.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/bary_entity_element.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/bary_grid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/bary_grid_view.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/bary_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/bisection.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/edge_creator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/edges_on_face.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/element_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/element_typeinfo.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/element_typelist.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/element_types.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entities.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entity_edge.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entity_element.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entity_face.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entity_iterator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entity_node.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/entity_pointer.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/generic_coordinate_transformator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/geometry_identity_mapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/geometry_mapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/geometry_sphere_mapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entities.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entity_edge.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entity_element.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entity_iterator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entity_node.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entity_pointer.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid2d_entity_seed.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid_creator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid_partitioner.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/grid_view.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/index_set.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_iterator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_base_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_mapping_data.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/intersection_set.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/make_entity_pair.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/mapper_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/no_intersection_filter.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/num_nodes_to_grid2d_element_type.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/orientation_container.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/orientation_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/permutation.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/surface_grid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/surface_grid_init_from_input.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/surface_grid_init_from_other.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/surface_grid_init_from_volume_grid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/surface_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/volume2d_grid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/volume2d_grid_init_from_input.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/volume2d_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/volume_to_surface_mapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/grid/volume_traits.hpp"
    )
endif()

