# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/excalibur/AdvNum/Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/excalibur/AdvNum/Code/build

# Include any dependencies generated for this target.
include HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/depend.make

# Include the progress variables for this target.
include HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/progress.make

# Include the compile flags for this target's objects.
include HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/flags.make

HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.o: HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/flags.make
HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.o: ../HMAT/HMAT_2D/is_admissible_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/excalibur/AdvNum/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.o"
	cd /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.o -c /home/excalibur/AdvNum/Code/HMAT/HMAT_2D/is_admissible_test.cpp

HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.i"
	cd /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/excalibur/AdvNum/Code/HMAT/HMAT_2D/is_admissible_test.cpp > CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.i

HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.s"
	cd /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/excalibur/AdvNum/Code/HMAT/HMAT_2D/is_admissible_test.cpp -o CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.s

# Object files for target is_admissible_test_solution
is_admissible_test_solution_OBJECTS = \
"CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.o"

# External object files for target is_admissible_test_solution
is_admissible_test_solution_EXTERNAL_OBJECTS =

bin/is_admissible_test_solution: HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/is_admissible_test.cpp.o
bin/is_admissible_test_solution: HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/build.make
bin/is_admissible_test_solution: HMAT/HMAT_2D/libhmat_2d_sols.so
bin/is_admissible_test_solution: HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/excalibur/AdvNum/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/is_admissible_test_solution"
	cd /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/is_admissible_test_solution.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/build: bin/is_admissible_test_solution

.PHONY : HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/build

HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/clean:
	cd /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D && $(CMAKE_COMMAND) -P CMakeFiles/is_admissible_test_solution.dir/cmake_clean.cmake
.PHONY : HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/clean

HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/depend:
	cd /home/excalibur/AdvNum/Code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/excalibur/AdvNum/Code /home/excalibur/AdvNum/Code/HMAT/HMAT_2D /home/excalibur/AdvNum/Code/build /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D /home/excalibur/AdvNum/Code/build/HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : HMAT/HMAT_2D/CMakeFiles/is_admissible_test_solution.dir/depend

