# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/zcy/Documents/projects/MarvelPhysics

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zcy/Documents/projects/MarvelPhysics/cmake

# Include any dependencies generated for this target.
include Point_Sys/src/CMakeFiles/PointSys.dir/depend.make

# Include the progress variables for this target.
include Point_Sys/src/CMakeFiles/PointSys.dir/progress.make

# Include the compile flags for this target's objects.
include Point_Sys/src/CMakeFiles/PointSys.dir/flags.make

Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o: Point_Sys/src/CMakeFiles/PointSys.dir/flags.make
Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o: ../Point_Sys/src/data_stream.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zcy/Documents/projects/MarvelPhysics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PointSys.dir/data_stream.cc.o -c /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/data_stream.cc

Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PointSys.dir/data_stream.cc.i"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/data_stream.cc > CMakeFiles/PointSys.dir/data_stream.cc.i

Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PointSys.dir/data_stream.cc.s"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/data_stream.cc -o CMakeFiles/PointSys.dir/data_stream.cc.s

Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.requires:

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.requires

Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.provides: Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.requires
	$(MAKE) -f Point_Sys/src/CMakeFiles/PointSys.dir/build.make Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.provides.build
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.provides

Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.provides.build: Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o


Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o: Point_Sys/src/CMakeFiles/PointSys.dir/flags.make
Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o: ../Point_Sys/src/gen_points.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zcy/Documents/projects/MarvelPhysics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PointSys.dir/gen_points.cc.o -c /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/gen_points.cc

Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PointSys.dir/gen_points.cc.i"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/gen_points.cc > CMakeFiles/PointSys.dir/gen_points.cc.i

Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PointSys.dir/gen_points.cc.s"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/gen_points.cc -o CMakeFiles/PointSys.dir/gen_points.cc.s

Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.requires:

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.requires

Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.provides: Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.requires
	$(MAKE) -f Point_Sys/src/CMakeFiles/PointSys.dir/build.make Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.provides.build
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.provides

Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.provides.build: Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o


Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o: Point_Sys/src/CMakeFiles/PointSys.dir/flags.make
Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o: ../Point_Sys/src/geometry.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zcy/Documents/projects/MarvelPhysics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PointSys.dir/geometry.cc.o -c /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/geometry.cc

Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PointSys.dir/geometry.cc.i"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/geometry.cc > CMakeFiles/PointSys.dir/geometry.cc.i

Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PointSys.dir/geometry.cc.s"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/geometry.cc -o CMakeFiles/PointSys.dir/geometry.cc.s

Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.requires:

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.requires

Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.provides: Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.requires
	$(MAKE) -f Point_Sys/src/CMakeFiles/PointSys.dir/build.make Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.provides.build
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.provides

Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.provides.build: Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o


Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o: Point_Sys/src/CMakeFiles/PointSys.dir/flags.make
Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o: ../Point_Sys/src/get_nn.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zcy/Documents/projects/MarvelPhysics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PointSys.dir/get_nn.cc.o -c /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/get_nn.cc

Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PointSys.dir/get_nn.cc.i"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/get_nn.cc > CMakeFiles/PointSys.dir/get_nn.cc.i

Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PointSys.dir/get_nn.cc.s"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/get_nn.cc -o CMakeFiles/PointSys.dir/get_nn.cc.s

Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.requires:

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.requires

Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.provides: Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.requires
	$(MAKE) -f Point_Sys/src/CMakeFiles/PointSys.dir/build.make Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.provides.build
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.provides

Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.provides.build: Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o


Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o: Point_Sys/src/CMakeFiles/PointSys.dir/flags.make
Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o: ../Point_Sys/src/points_energy.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zcy/Documents/projects/MarvelPhysics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PointSys.dir/points_energy.cc.o -c /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/points_energy.cc

Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PointSys.dir/points_energy.cc.i"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/points_energy.cc > CMakeFiles/PointSys.dir/points_energy.cc.i

Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PointSys.dir/points_energy.cc.s"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src/points_energy.cc -o CMakeFiles/PointSys.dir/points_energy.cc.s

Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.requires:

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.requires

Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.provides: Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.requires
	$(MAKE) -f Point_Sys/src/CMakeFiles/PointSys.dir/build.make Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.provides.build
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.provides

Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.provides.build: Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o


# Object files for target PointSys
PointSys_OBJECTS = \
"CMakeFiles/PointSys.dir/data_stream.cc.o" \
"CMakeFiles/PointSys.dir/gen_points.cc.o" \
"CMakeFiles/PointSys.dir/geometry.cc.o" \
"CMakeFiles/PointSys.dir/get_nn.cc.o" \
"CMakeFiles/PointSys.dir/points_energy.cc.o"

# External object files for target PointSys
PointSys_EXTERNAL_OBJECTS =

lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o
lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o
lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o
lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o
lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o
lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/build.make
lib/libPointSys.so: /usr/lib/x86_64-linux-gnu/libtbb.so
lib/libPointSys.so: Point_Sys/src/CMakeFiles/PointSys.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zcy/Documents/projects/MarvelPhysics/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX shared library ../../lib/libPointSys.so"
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PointSys.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Point_Sys/src/CMakeFiles/PointSys.dir/build: lib/libPointSys.so

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/build

Point_Sys/src/CMakeFiles/PointSys.dir/requires: Point_Sys/src/CMakeFiles/PointSys.dir/data_stream.cc.o.requires
Point_Sys/src/CMakeFiles/PointSys.dir/requires: Point_Sys/src/CMakeFiles/PointSys.dir/gen_points.cc.o.requires
Point_Sys/src/CMakeFiles/PointSys.dir/requires: Point_Sys/src/CMakeFiles/PointSys.dir/geometry.cc.o.requires
Point_Sys/src/CMakeFiles/PointSys.dir/requires: Point_Sys/src/CMakeFiles/PointSys.dir/get_nn.cc.o.requires
Point_Sys/src/CMakeFiles/PointSys.dir/requires: Point_Sys/src/CMakeFiles/PointSys.dir/points_energy.cc.o.requires

.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/requires

Point_Sys/src/CMakeFiles/PointSys.dir/clean:
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src && $(CMAKE_COMMAND) -P CMakeFiles/PointSys.dir/cmake_clean.cmake
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/clean

Point_Sys/src/CMakeFiles/PointSys.dir/depend:
	cd /home/zcy/Documents/projects/MarvelPhysics/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zcy/Documents/projects/MarvelPhysics /home/zcy/Documents/projects/MarvelPhysics/Point_Sys/src /home/zcy/Documents/projects/MarvelPhysics/cmake /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src /home/zcy/Documents/projects/MarvelPhysics/cmake/Point_Sys/src/CMakeFiles/PointSys.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Point_Sys/src/CMakeFiles/PointSys.dir/depend

