include(ExternalProject) # To download Eigen if it is missing

find_package(EIGEN3)

if(${EIGEN3_FOUND})

    include_directories(${EIGEN3_INCLUDE_DIR})

else()

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

