# Install script for directory: /home/excalibur/AdvNum/Code/third_party/Betl2/Library

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/cmdl_parser/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/input_interface/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/gmsh_input/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/linalg/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/grid/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/fe/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/utils/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/geometry/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/quadrature/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_integration/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/integration/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/fundsol/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/functional/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/analytical_functions/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/vtu_exporter/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/math_special_functions/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/material/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_operator/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/sparse_operators/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/fem_operator/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_symmetries/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/multilevel/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/hypre/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/pims/cmake_install.cmake")
  include("/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/driver/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
