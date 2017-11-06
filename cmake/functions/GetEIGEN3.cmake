include(ExternalProject) # To download Eigen if it is missing

find_package(EIGEN3)

# Added the following requirement for BEtl2 to work
#if(EIGEN3_VERSION VERSION_GREATER 3.2.7 OR EIGEN3_VERSION VERSION_EQUAL 3.2.7)
#  message("Eigen version >= 3.2.7 is not supported. Downloading 3.2.7 instead.")
#  set(EIGEN3_FOUND FALSE)
#endif()
if(EIGEN3_VERSION VERSION_LESS 3.3.0)
  message("Old Eigen version. Downloading 3.3.4 instead.")
  set(EIGEN3_FOUND FALSE)
endif()

if(${EIGEN3_FOUND})

  include_directories(${EIGEN3_INCLUDE_DIR})
  add_custom_target(Eigen) # dependency dummy

else()
  #  if not found system wide download. Changed eigen version from 3.3.4 to 3.2.7
  # for Betl2 to work
    SET(DOWNLOADING_EIGEN ON)
    message("-- Downloading Eigen3")
    ExternalProject_Add(
        Eigen
        URL http://bitbucket.org/eigen/eigen/get/3.3.4.zip
         SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
        DOWNLOAD_NO_PROGRESS 1
        CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install)
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)

endif()

