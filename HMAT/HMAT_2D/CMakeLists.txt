cmake_minimum_required(VERSION 3.1)

project(HMAT_2D)

add_subdirectory(include)
add_subdirectory(src)

# Main target to get familiar with HMAT_2D code
# Template target
add_library(hmat_2d_tmps SHARED ${headers} ${cpp_files} ${cpp_files_tmps})
add_executable(low_rank_test_2D_template low_rank_test.cpp)
target_link_libraries(low_rank_test_2D_template hmat_2d_tmps)

# Solution target
add_library(hmat_2d_sols SHARED ${headers} ${cpp_files} ${cpp_files_sols})
add_executable(low_rank_test_2D_solution low_rank_test.cpp)
target_link_libraries(low_rank_test_2D_solution hmat_2d_sols)

# Test target for cluster tree
# Template target
add_executable(ctree_test_template ctree_test.cpp)
target_link_libraries(ctree_test_template hmat_2d_tmps)

# Solution target
add_executable(ctree_test_solution ctree_test.cpp)
target_link_libraries(ctree_test_solution hmat_2d_sols)

# Test target for admissibility condition
# Template target
add_executable(is_admissible_test_template is_admissible_test.cpp)
target_link_libraries(is_admissible_test_template hmat_2d_tmps)

# Solution target
add_executable(is_admissible_test_solution is_admissible_test.cpp)
target_link_libraries(is_admissible_test_solution hmat_2d_sols)

# Test target for separation into far- and near-field
# Template target
add_executable(nearfar_test_template nearfar_test.cpp)
target_link_libraries(nearfar_test_template hmat_2d_tmps)

# Solution target
add_executable(nearfar_test_solution nearfar_test.cpp)
target_link_libraries(nearfar_test_solution hmat_2d_sols)

# Test target for V-matrix per node
# Template target
add_executable(vnode_matrix_test_template vnode_matrix_test.cpp)
target_link_libraries(vnode_matrix_test_template hmat_2d_tmps)

# Solution target
add_executable(vnode_matrix_test_solution vnode_matrix_test.cpp)
target_link_libraries(vnode_matrix_test_solution hmat_2d_sols)

