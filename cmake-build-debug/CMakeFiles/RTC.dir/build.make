# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = E:\apps\CLion\ch-0\203.7148.70\bin\cmake\win\bin\cmake.exe

# The command to remove a file.
RM = E:\apps\CLion\ch-0\203.7148.70\bin\cmake\win\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = F:\RTC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = F:\RTC\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/RTC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/RTC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/RTC.dir/flags.make

CMakeFiles/RTC.dir/main.cpp.obj: CMakeFiles/RTC.dir/flags.make
CMakeFiles/RTC.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=F:\RTC\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/RTC.dir/main.cpp.obj"
	E:\CodeBlocks\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\RTC.dir\main.cpp.obj -c F:\RTC\main.cpp

CMakeFiles/RTC.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RTC.dir/main.cpp.i"
	E:\CodeBlocks\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E F:\RTC\main.cpp > CMakeFiles\RTC.dir\main.cpp.i

CMakeFiles/RTC.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RTC.dir/main.cpp.s"
	E:\CodeBlocks\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S F:\RTC\main.cpp -o CMakeFiles\RTC.dir\main.cpp.s

# Object files for target RTC
RTC_OBJECTS = \
"CMakeFiles/RTC.dir/main.cpp.obj"

# External object files for target RTC
RTC_EXTERNAL_OBJECTS =

RTC.exe: CMakeFiles/RTC.dir/main.cpp.obj
RTC.exe: CMakeFiles/RTC.dir/build.make
RTC.exe: CMakeFiles/RTC.dir/linklibs.rsp
RTC.exe: CMakeFiles/RTC.dir/objects1.rsp
RTC.exe: CMakeFiles/RTC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=F:\RTC\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable RTC.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\RTC.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/RTC.dir/build: RTC.exe

.PHONY : CMakeFiles/RTC.dir/build

CMakeFiles/RTC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\RTC.dir\cmake_clean.cmake
.PHONY : CMakeFiles/RTC.dir/clean

CMakeFiles/RTC.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" F:\RTC F:\RTC F:\RTC\cmake-build-debug F:\RTC\cmake-build-debug F:\RTC\cmake-build-debug\CMakeFiles\RTC.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/RTC.dir/depend

