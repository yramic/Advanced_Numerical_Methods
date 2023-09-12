# Provides variable PROBLEM_NAME
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_variables.cmake)

# Provides functions build_problem and build_test
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_rules.cmake)

# pass correct arguemnts to build rules
function(build PROBLEM_NAME DIR SOLUTION)
    set(PROBLEM_TARGET ${PROBLEM_NAME}_${DIR})
    set(TEST_TARGET ${PROBLEM_NAME}_test_${DIR})

    # problem
    build_problem(${PROBLEM_TARGET} ${DIR} ${PROBLEM_TARGET})
    target_compile_definitions(${PROBLEM_TARGET} PRIVATE SOLUTION=${SOLUTION})
    target_compile_definitions(${PROBLEM_TARGET}.static PRIVATE SOLUTION=${SOLUTION})

    # tests
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/test)
        build_test(${TEST_TARGET} ${PROBLEM_TARGET} ${DIR} ${TEST_TARGET})
        target_compile_definitions(${TEST_TARGET} PRIVATE SOLUTION=${SOLUTION})
    endif()
endfunction(build)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/all)
    build(${PROBLEM_NAME} all 1)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mastersolution)
    build(${PROBLEM_NAME} mastersolution 1)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mysolution)
    build(${PROBLEM_NAME} mysolution 0)
endif()
