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
include BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/depend.make

# Include the progress variables for this target.
include BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/progress.make

# Include the compile flags for this target's objects.
include BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/flags.make

BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/test_build.cpp.o: BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/flags.make
BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/test_build.cpp.o: ../BEM/CppHilbert/Examples/test_build.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/excalibur/AdvNum/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/test_build.cpp.o"
	cd /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cmtest.dir/test_build.cpp.o -c /home/excalibur/AdvNum/Code/BEM/CppHilbert/Examples/test_build.cpp

BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/test_build.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cmtest.dir/test_build.cpp.i"
	cd /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/excalibur/AdvNum/Code/BEM/CppHilbert/Examples/test_build.cpp > CMakeFiles/cmtest.dir/test_build.cpp.i

BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/test_build.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cmtest.dir/test_build.cpp.s"
	cd /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/excalibur/AdvNum/Code/BEM/CppHilbert/Examples/test_build.cpp -o CMakeFiles/cmtest.dir/test_build.cpp.s

# Object files for target cmtest
cmtest_OBJECTS = \
"CMakeFiles/cmtest.dir/test_build.cpp.o"

# External object files for target cmtest
cmtest_EXTERNAL_OBJECTS =

bin/cmtest: BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/test_build.cpp.o
bin/cmtest: BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/build.make
bin/cmtest: BEM/CppHilbert/Library/libCppHilbert.so
bin/cmtest: BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/excalibur/AdvNum/Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../bin/cmtest"
	cd /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/build: bin/cmtest

.PHONY : BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/build

BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/clean:
	cd /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples && $(CMAKE_COMMAND) -P CMakeFiles/cmtest.dir/cmake_clean.cmake
.PHONY : BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/clean

BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/depend:
	cd /home/excalibur/AdvNum/Code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/excalibur/AdvNum/Code /home/excalibur/AdvNum/Code/BEM/CppHilbert/Examples /home/excalibur/AdvNum/Code/build /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples /home/excalibur/AdvNum/Code/build/BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : BEM/CppHilbert/Examples/CMakeFiles/cmtest.dir/depend

