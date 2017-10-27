cmake_minimum_required(VERSION 3.1)

project(HMAT_2D)

add_subdirectory(include)
add_subdirectory(src)

add_executable(low_rank_test
               low_rank_test.cpp ${cpp_files} ${headers})

add_executable(ctree_test
               ctree_test.cpp ${cpp_files} ${headers})

add_executable(vnode_matrix_test
               vnode_matrix_test.cpp ${cpp_files} ${headers})

add_executable(is_admissible_test
               is_admissible_test.cpp ${cpp_files} ${headers})

add_executable(nearfar_test
               nearfar_test.cpp ${cpp_files} ${headers})
