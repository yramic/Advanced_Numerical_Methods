cmake_minimum_required( VERSION 2.8 )

project(parametricbem)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(DIRECTORIES src test)
#link_directories(${PROJECT_BINARY_DIR}/src)

foreach(directory ${DIRECTORIES})
	add_subdirectory(${directory})
endforeach(directory)

# This needs to be in both the top-level and the test CMake file.
enable_testing()

# No main target to build..
