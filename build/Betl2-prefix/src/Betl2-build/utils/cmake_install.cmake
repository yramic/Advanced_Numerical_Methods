# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libutils.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/utils/libutils.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/collision_detection.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/deleter.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/detail_matio_exporter.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/detail_matio_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/element_metrics.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/element_metrics_entity_distance.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/element_metrics_global_mesh_size.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/element_metrics_local_mesh_size.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/element_metrics_mesh_regularity.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/if_then_else.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/make_matrix.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/matiostream.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/matrix_eigen_allocator.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/matrix_operations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/matrix_storage_scheme.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/num_threads.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/numeric_traits.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/set_zero.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/sparsestream.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/utils/static_iterate.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/utils" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/collision_detection.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/deleter.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/detail_matio_exporter.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/detail_matio_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/element_metrics.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/element_metrics_entity_distance.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/element_metrics_global_mesh_size.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/element_metrics_local_mesh_size.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/element_metrics_mesh_regularity.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/if_then_else.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/make_matrix.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/matiostream.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/matrix_eigen_allocator.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/matrix_operations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/matrix_storage_scheme.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/num_threads.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/numeric_traits.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/set_zero.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/sparsestream.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/utils/static_iterate.hpp"
    )
endif()

