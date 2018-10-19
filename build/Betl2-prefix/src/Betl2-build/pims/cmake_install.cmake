# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims

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
   "/home/excalibur/AdvNum/Code/build/betl2_install/lib/libpetsc_bindings.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/lib" TYPE STATIC_LIBRARY FILES "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/pims/libpetsc_bindings.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/converter.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/combined_preconditioner.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/error.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/export_matrix.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/initialize.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/matrix_wrapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/matrix_shells.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/monitoring_functions.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/operations.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/preconditioner.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/vector_wrapper.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/world.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/solver.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/bem_preconditioner.hpp;/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc/sparsity_policy_petsc.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc" TYPE FILE FILES
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/converter.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/combined_preconditioner.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/error.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/export_matrix.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/initialize.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/matrix_wrapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/matrix_shells.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/monitoring_functions.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/operations.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/preconditioner.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/vector_wrapper.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/world.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/solver.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/bem_preconditioner.hpp"
    "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc/sparsity_policy_petsc.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/excalibur/AdvNum/Code/build/betl2_install/include/pims/petsc.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/excalibur/AdvNum/Code/build/betl2_install/include/pims" TYPE FILE FILES "/home/excalibur/AdvNum/Code/third_party/Betl2/Library/pims/petsc.hpp")
endif()

