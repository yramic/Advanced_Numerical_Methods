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
CMAKE_SOURCE_DIR = /home/excalibur/AdvNum/Code/third_party/Betl2/Library

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build

# Include any dependencies generated for this target.
include bem_operator/CMakeFiles/bem_operator.dir/depend.make

# Include the progress variables for this target.
include bem_operator/CMakeFiles/bem_operator.dir/progress.make

# Include the compile flags for this target's objects.
include bem_operator/CMakeFiles/bem_operator.dir/flags.make

# Object files for target bem_operator
bem_operator_OBJECTS =

# External object files for target bem_operator
bem_operator_EXTERNAL_OBJECTS =

bem_operator/libbem_operator.a: bem_operator/CMakeFiles/bem_operator.dir/build.make
bem_operator/libbem_operator.a: bem_operator/CMakeFiles/bem_operator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX static library libbem_operator.a"
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_operator && $(CMAKE_COMMAND) -P CMakeFiles/bem_operator.dir/cmake_clean_target.cmake
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_operator && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bem_operator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
bem_operator/CMakeFiles/bem_operator.dir/build: bem_operator/libbem_operator.a

.PHONY : bem_operator/CMakeFiles/bem_operator.dir/build

bem_operator/CMakeFiles/bem_operator.dir/clean:
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_operator && $(CMAKE_COMMAND) -P CMakeFiles/bem_operator.dir/cmake_clean.cmake
.PHONY : bem_operator/CMakeFiles/bem_operator.dir/clean

bem_operator/CMakeFiles/bem_operator.dir/depend:
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/excalibur/AdvNum/Code/third_party/Betl2/Library /home/excalibur/AdvNum/Code/third_party/Betl2/Library/bem_operator /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_operator /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/bem_operator/CMakeFiles/bem_operator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : bem_operator/CMakeFiles/bem_operator.dir/depend

