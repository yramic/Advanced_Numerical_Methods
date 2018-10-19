# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/MACRO_make_quadrature.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/MACRO_make_quadrature_segment_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_data.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_hexa.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_hexa_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_list.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_prism.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_prism_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_pyramid.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_pyramid_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_quad.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_quad_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_segment.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_segment_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_tetra.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_tetra_impl.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_tria.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/quadrature_tria_impl.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/MACRO_make_quadrature.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/MACRO_make_quadrature_segment_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_data.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_hexa.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_hexa_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_list.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_prism.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_prism_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_pyramid.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_pyramid_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_quad.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_quad_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_segment.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_segment_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_tetra.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_tetra_impl.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_tria.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/quadrature_tria_impl.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature/segment_data")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/quadrature" TYPE DIRECTORY FILES "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/quadrature/segment_data")
endif()

