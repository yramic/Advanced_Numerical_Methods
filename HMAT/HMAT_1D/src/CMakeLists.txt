set(cpp_files
    ${CMAKE_CURRENT_LIST_DIR}/block_cluster.cpp
    ${CMAKE_CURRENT_LIST_DIR}/block_nearf.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cheby.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ctree.cpp
    ${CMAKE_CURRENT_LIST_DIR}/hierarchical_partition.cpp
    ${CMAKE_CURRENT_LIST_DIR}/is_admissible.cpp
    ${CMAKE_CURRENT_LIST_DIR}/kernel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/low_rank_app.cpp
    ${CMAKE_CURRENT_LIST_DIR}/node.cpp
    ${CMAKE_CURRENT_LIST_DIR}/point.cpp
    PARENT_SCOPE
)

add_subdirectory(solutions)
set(  cpp_files_sols
    ${cpp_files_sols}
    PARENT_SCOPE
)
add_subdirectory(templates)
set(  cpp_files_tmps
    ${cpp_files_tmps}
    PARENT_SCOPE
)

add_subdirectory(uni-direct)
set(  cpp_files_uni
    ${cpp_files_uni}
    PARENT_SCOPE
)
set(  cpp_files_uni_sols
    ${cpp_files_uni_sols}
    PARENT_SCOPE
)
set(  cpp_files_uni_tmps
    ${cpp_files_uni_tmps}
    PARENT_SCOPE
)

