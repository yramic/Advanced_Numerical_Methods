add_executable(clustertree_revised clustertree.h clustertree.cpp clustertree_main.cpp)
add_executable(matrixpartition_revised clustertree.h clustertree.cpp matrixpartition.h
  matrixpartition.cpp matrixpartition_main.cpp)

target_link_libraries(clustertree_revised PUBLIC Eigen3::Eigen)
target_link_libraries(matrixpartition_revised PUBLIC Eigen3::Eigen)

add_library(Clustering_revised clustertree.cpp matrixpartition.cpp)
target_link_libraries(Clustering_revised Eigen3::Eigen)
