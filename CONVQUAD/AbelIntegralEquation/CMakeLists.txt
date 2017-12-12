add_executable_advnumcse(solution AbelIntegralEquation.cpp)
get_target_name_advnumcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_advnumcse(template AbelIntegralEquation.cpp)
get_target_name_advnumcse(template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
