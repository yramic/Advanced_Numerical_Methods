cmake_minimum_required( VERSION 2.8 )

set(CMAKE_CXX_FLAGS "-std=c++11")
set(DIRECTORIES src)

foreach(directory ${DIRECTORIES})
	add_subdirectory(${directory})
endforeach(directory)

