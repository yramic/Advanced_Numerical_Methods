# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libhypre_bindings.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/hypre/libhypre_bindings.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/hypre.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/hypre_error.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/hypre_initialize.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/hypre_matrix_wrapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/hypre_vector_wrapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/hypre_world.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre/sparsity_policy_hypre.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/hypre" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/hypre.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/hypre_error.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/hypre_initialize.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/hypre_matrix_wrapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/hypre_vector_wrapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/hypre_world.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/hypre/sparsity_policy_hypre.hpp"
    )
endif()

