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
include driver/CMakeFiles/driver.dir/depend.make

# Include the progress variables for this target.
include driver/CMakeFiles/driver.dir/progress.make

# Include the compile flags for this target's objects.
include driver/CMakeFiles/driver.dir/flags.make

# Object files for target driver
driver_OBJECTS =

# External object files for target driver
driver_EXTERNAL_OBJECTS =

driver/libdriver.a: driver/CMakeFiles/driver.dir/build.make
driver/libdriver.a: driver/CMakeFiles/driver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX static library libdriver.a"
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/driver && $(CMAKE_COMMAND) -P CMakeFiles/driver.dir/cmake_clean_target.cmake
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/driver && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/driver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
driver/CMakeFiles/driver.dir/build: driver/libdriver.a

.PHONY : driver/CMakeFiles/driver.dir/build

driver/CMakeFiles/driver.dir/clean:
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/driver && $(CMAKE_COMMAND) -P CMakeFiles/driver.dir/cmake_clean.cmake
.PHONY : driver/CMakeFiles/driver.dir/clean

driver/CMakeFiles/driver.dir/depend:
	cd /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/excalibur/AdvNum/Code/third_party/Betl2/Library /home/excalibur/AdvNum/Code/third_party/Betl2/Library/driver /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/driver /home/excalibur/AdvNum/Code/build/Betl2-prefix/src/Betl2-build/driver/CMakeFiles/driver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : driver/CMakeFiles/driver.dir/depend

