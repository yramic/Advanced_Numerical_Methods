add_executable(clustertree clustertree.h clustertree.cpp clustertree_main.cpp)
add_executable(matrixpartition clustertree.h clustertree.cpp matrixpartition.h
  matrixpartition.cpp matrixpartition_main.cpp)

target_link_libraries(clustertree PUBLIC Eigen3::Eigen)
target_link_libraries(matrixpartition PUBLIC Eigen3::Eigen)

add_library(Clustering clustertree.cpp matrixpartition.cpp)
target_link_libraries(Clustering Eigen3::Eigen)
