add_subdirectory(include)
add_subdirectory(src)

add_executable(low_rank_test low_rank_test.cpp ${cpp_files} ${headers})

