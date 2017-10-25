add_executable_advnumcse(solution DBEMmain.cpp)
get_target_name_advnumcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)
target_link_libraries(${SOLUTION} CppHilbert)

add_executable_advnumcse(template DBEMmain.cpp)
get_target_name_advnumcse(template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
target_link_libraries(${TEMPLATE} CppHilbert)
