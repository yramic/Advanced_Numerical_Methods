# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_collection.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection_hex.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/geometry_impl_intersection_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_segment.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gradient_shape_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_precalc.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_segment.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/shape_functions_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/integration_element.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry/gram_and_normal.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/geometry" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_collection.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection_hex.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/geometry_impl_intersection_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_segment.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gradient_shape_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_precalc.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_segment.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/shape_functions_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/integration_element.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/geometry/gram_and_normal.hpp"
    )
endif()

