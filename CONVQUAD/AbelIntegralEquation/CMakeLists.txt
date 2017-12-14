# Template target
add_executable(convquad_abel-integral-equation_template AbelIntegralEquation.cpp)
target_link_libraries(convquad_abel-integral-equation_template CppHilbert)
target_compile_definitions(convquad_abel-integral-equation_template PRIVATE SOLUTION=0 INTERNAL=0)

# Solution target
add_executable(convquad_abel-integral-equation_solution AbelIntegralEquation.cpp)
target_link_libraries(convquad_abel-integral-equation_solution CppHilbert)
target_compile_definitions(convquad_abel-integral-equation_solution PRIVATE SOLUTION=1 INTERNAL=1)

