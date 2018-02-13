include(CMake/functions/get_target_name_advnumcse.cmake)

macro(LIST_REPLACE LIST INDEX NEWVALUE)
    list(INSERT ${LIST} ${INDEX} ${NEWVALUE})
    MATH(EXPR __INDEX "${INDEX} + 1")
    list(REMOVE_AT ${LIST} ${__INDEX})
endmacro(LIST_REPLACE)

# note: install rpath is disabled for now since this is only required for relocatable binaries

macro(add_executable_advnumcse_betl2 _name)
    # cmake targets need to be unique
    # as such we prepend a prefix from the source path
    # all of this is done by get_executable_name_advnumcse
    get_target_name_advnumcse(${_name} target_name)
    # reassemble argv for add_executable macro with new target name
    set(argv_parsed ${ARGV})
    LIST_REPLACE(argv_parsed 0 "${target_name}")

    # create directory that will contain the binary
    set(TARGET_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # retrieve relative path of runtime output directory (e.g., the directory to which the target is compiled to)
    string(REPLACE "${CMAKE_SOURCE_DIR}" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" TARGET_RUNTIME_OUTPUT_DIRECTORY ${TARGET_RUNTIME_OUTPUT_DIRECTORY})

    # get relative path to binary dir (used as an rpath)
#    string(REPLACE "${CMAKE_BINARY_DIR}" "" RELATIVE_PATH_TO_BINARY_DIR "${TARGET_RUNTIME_OUTPUT_DIRECTORY}")
#    string(REGEX REPLACE "/[^/]+" "../" RELATIVE_PATH_TO_BINARY_DIR "${RELATIVE_PATH_TO_BINARY_DIR}")

    # invoke built-in add_executable
    add_executable(${argv_parsed})
    # link with mathgl and figure class
    if(TARGET ${target_name})
        set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}"
        OUTPUT_NAME "${_name}"
#        INSTALL_RPATH "@executable_path/${RELATIVE_PATH_TO_BINARY_DIR}/mathgl_install/lib"
)
	target_link_libraries(${target_name} ${BETL2_LIBRARIES})
       	add_dependencies(${target_name} Eigen)
        add_dependencies(${target_name} Betl2)
    endif()
endmacro()

