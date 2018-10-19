################################################################################
#
# Find ETH headers and lib
#
# ETH_INCLUDEDIR	 - prefix for header files
# ETH_LIBRARYDIR    - library path
#
# ETH_FOUND       - true when found
################################################################################

#check if ETH_ROOT is defined in envir and use it
if (NOT ETH_ROOT AND NOT $ENV{ETH_ROOT} STREQUAL "")
  set (ETH_ROOT $ENV{ETH_ROOT})
endif ()

# If ETH_INCLUDEDIR was defined in the environment, use it.
IF( NOT $ENV{ETH_INCLUDEDIR} STREQUAL "" )
  set(ETH_INCLUDEDIR $ENV{ETH_INCLUDEDIR})
ENDIF( NOT $ENV{ETH_INCLUDEDIR} STREQUAL "" )

# If ETH_LIBRARYDIR was defined in the environment, use it.
IF( NOT $ENV{ETH_LIBRARYDIR} STREQUAL "" )
  set(ETH_LIBRARYDIR $ENV{ETH_LIBRARYDIR})
ENDIF( NOT $ENV{ETH_LIBRARYDIR} STREQUAL "" )

# figure out what the next statement does
IF( ETH_ROOT )
   file(TO_CMAKE_PATH ${ETH_ROOT} ETH_ROOT)
ENDIF( ETH_ROOT )

FIND_PATH(ETH_INCLUDE_DIRS cmdl_parser/cmdl_parser.hpp 
  ${ETH_ROOT}
  ${ETH_INCLUDEDIR}
  ${ETH_ROOT}/include
  ${ETH_ROOT}/Include
  /usr/local/include
  "$ENV{ProgramFiles}/include"
  "$ENV{ProgramW6432}/include"
)

SET (ETH_LIBRARY_DIRS
  ${ETH_ROOT}
  ${ETH_ROOT}/Lib
  ${ETH_ROOT}/lib
  /usr/local/lib
  "$ENV{ProgramFiles}/lib"
  "$ENV{ProgramW6432}/lib"
)

SET( LIBNAMES
  cmdl_parser
  gmsh_input
  eth_base
  interface
  
  )

FIND_LIBRARY(ETH_BASE_LIB           eth_base       ${ETH_LIBRARY_DIRS})
FIND_LIBRARY(ETH_CMDL_PARSER_LIB    cmdl_parser    ${ETH_LIBRARY_DIRS})
FIND_LIBRARY(ETH_GMSH_INPUT_LIB     gmsh_input     ${ETH_LIBRARY_DIRS})
FIND_LIBRARY(ETH_BASE_LIB			base		   ${ETH_LIBRARY_DIRS})
FIND_LIBRARY(ETH_INTERFACE_LIB		interface	   ${ETH_LIBRARY_DIRS})

message (STATUS ${ETH_INCLUDE_DIRS})
message (STATUS ${ETH_BASE_LIB})
message (STATUS ${ETH_CMDL_PARSER_LIB})
message (STATUS ${ETH_GMSH_INPUT_LIB})
MESSAGE (STATUS ${ETH_BASE_LIB})
MESSAGE (STATUS ${ETH_INTERFACE_LIB})

IF(ETH_INCLUDE_DIRS AND ETH_BASE_LIB AND ETH_CMDL_PARSER_LIB AND ETH_GMSH_INPUT_LIB )
  SET( ETH_FOUND "YES")
ENDIF()

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  ETH
  DEFAULT_MSG
  ETH_INCLUDE_DIRS
  ETH_BASE_LIB
  ETH_CMDL_PARSER_LIB
  ETH_GMSH_INPUT_LIB
  ETH_BASE_LIB
  ETH_INTERFACE_LIB
  )

MARK_AS_ADVANCED( ETH_INCLUDE_DIRS
  ETH_BASE_LIB
  ETH_CMDL_PARSER_LIB
  ETH_GMSH_INPUT_LIB
  )

SET( ETH_LIBRARIES
  ${ETH_BASE_LIB}
  ${ETH_CMDL_PARSER_LIB}
  ${ETH_GMSH_INPUT_LIB}
  ${ETH_BASE_LIB}
  ${ETH_INTERFACE_LIB}
  )


