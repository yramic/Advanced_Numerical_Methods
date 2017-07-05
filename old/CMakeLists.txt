set(headers
    ${CMAKE_CURRENT_LIST_DIR}/cluster.hpp
)

set(cpp_files
    ${CMAKE_CURRENT_LIST_DIR}/changeable.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cluster.cpp
)

add_executable(low_rank_old main.cpp ${cpp_files} ${headers})

