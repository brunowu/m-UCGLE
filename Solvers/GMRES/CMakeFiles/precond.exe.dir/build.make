# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.6.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.6.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/xwu/phd_xwu/UCMGEL

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/xwu/phd_xwu/UCMGEL

# Include any dependencies generated for this target.
include Solvers/GMRES/CMakeFiles/precond.exe.dir/depend.make

# Include the progress variables for this target.
include Solvers/GMRES/CMakeFiles/precond.exe.dir/progress.make

# Include the compile flags for this target's objects.
include Solvers/GMRES/CMakeFiles/precond.exe.dir/flags.make

Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o: Solvers/GMRES/CMakeFiles/precond.exe.dir/flags.make
Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o: Solvers/GMRES/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xwu/phd_xwu/UCMGEL/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o"
	cd /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES && /Users/xwu/software/openmp/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/precond.exe.dir/main.cpp.o -c /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES/main.cpp

Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/precond.exe.dir/main.cpp.i"
	cd /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES && /Users/xwu/software/openmp/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES/main.cpp > CMakeFiles/precond.exe.dir/main.cpp.i

Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/precond.exe.dir/main.cpp.s"
	cd /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES && /Users/xwu/software/openmp/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES/main.cpp -o CMakeFiles/precond.exe.dir/main.cpp.s

Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.requires:

.PHONY : Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.requires

Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.provides: Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.requires
	$(MAKE) -f Solvers/GMRES/CMakeFiles/precond.exe.dir/build.make Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.provides.build
.PHONY : Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.provides

Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.provides.build: Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o


# Object files for target precond.exe
precond_exe_OBJECTS = \
"CMakeFiles/precond.exe.dir/main.cpp.o"

# External object files for target precond.exe
precond_exe_EXTERNAL_OBJECTS =

Solvers/GMRES/precond.exe: Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o
Solvers/GMRES/precond.exe: Solvers/GMRES/CMakeFiles/precond.exe.dir/build.make
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraext.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetrainout.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkostsqr.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraclassiclinalg.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraclassicnodeapi.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraclassic.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libepetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoskernels.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoskokkoscomm.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoskokkoscompat.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosremainder.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosnumerics.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoscomm.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosparameterlist.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoscore.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkosalgorithms.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoscontainers.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoscore.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoskokkoscomm.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoskokkoscompat.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosremainder.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosnumerics.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoscomm.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosparameterlist.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoscore.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoscore.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libbelostpetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libbelosepetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libbelos.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libaztecoo.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libthyratpetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libthyraepetraext.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libthyraepetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libthyracore.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libepetraext.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraext.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetrainout.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkostsqr.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraclassiclinalg.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraclassicnodeapi.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtpetraclassic.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libtriutils.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libepetra.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/librtop.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoskernels.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoskokkoscomm.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoskokkoscompat.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosremainder.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosnumerics.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoscomm.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchosparameterlist.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libteuchoscore.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkosalgorithms.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoscontainers.a
Solvers/GMRES/precond.exe: /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2/lib/libkokkoscore.a
Solvers/GMRES/precond.exe: /usr/lib/liblapack.dylib
Solvers/GMRES/precond.exe: /usr/lib/libblas.dylib
Solvers/GMRES/precond.exe: /usr/lib/libdl.dylib
Solvers/GMRES/precond.exe: Solvers/GMRES/CMakeFiles/precond.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/xwu/phd_xwu/UCMGEL/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable precond.exe"
	cd /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/precond.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Solvers/GMRES/CMakeFiles/precond.exe.dir/build: Solvers/GMRES/precond.exe

.PHONY : Solvers/GMRES/CMakeFiles/precond.exe.dir/build

Solvers/GMRES/CMakeFiles/precond.exe.dir/requires: Solvers/GMRES/CMakeFiles/precond.exe.dir/main.cpp.o.requires

.PHONY : Solvers/GMRES/CMakeFiles/precond.exe.dir/requires

Solvers/GMRES/CMakeFiles/precond.exe.dir/clean:
	cd /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES && $(CMAKE_COMMAND) -P CMakeFiles/precond.exe.dir/cmake_clean.cmake
.PHONY : Solvers/GMRES/CMakeFiles/precond.exe.dir/clean

Solvers/GMRES/CMakeFiles/precond.exe.dir/depend:
	cd /Users/xwu/phd_xwu/UCMGEL && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xwu/phd_xwu/UCMGEL /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES /Users/xwu/phd_xwu/UCMGEL /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES /Users/xwu/phd_xwu/UCMGEL/Solvers/GMRES/CMakeFiles/precond.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Solvers/GMRES/CMakeFiles/precond.exe.dir/depend
