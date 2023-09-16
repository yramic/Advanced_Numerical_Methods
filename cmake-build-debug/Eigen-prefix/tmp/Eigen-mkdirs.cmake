# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen"
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/src/Eigen-build"
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen_install"
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/tmp"
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/src/Eigen-stamp"
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/src"
  "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/src/Eigen-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/src/Eigen-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/Eigen-prefix/src/Eigen-stamp${cfgdir}") # cfgdir has leading slash
endif()
