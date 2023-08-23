add_executable_advnumcse(clustertree clustertree.h clustertree.cpp clustertree_main.cpp)
add_executable_advnumcse(matrixpartition clustertree.h clustertree.cpp matrixpartition.h
  matrixpartition.cpp matrixpartition_main.cpp)
add_executable_advnumcse(locallowrank clustertree.h clustertree.cpp matrixpartition.h
  matrixpartition.cpp locallowrank.h locallowrank.cpp locallowrank_main.cpp)
add_executable_advnumcse(demo Clustering.cpp clustertree.h clustertree.cpp
  matrixpartition.h matrixpartition.cpp)
