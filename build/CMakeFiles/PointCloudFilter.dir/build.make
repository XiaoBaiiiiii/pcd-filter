# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/vboxuser/code/pointcloudfilter

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vboxuser/code/pointcloudfilter/build

# Include any dependencies generated for this target.
include CMakeFiles/PointCloudFilter.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PointCloudFilter.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PointCloudFilter.dir/flags.make

CMakeFiles/PointCloudFilter.dir/main.cc.o: CMakeFiles/PointCloudFilter.dir/flags.make
CMakeFiles/PointCloudFilter.dir/main.cc.o: ../main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vboxuser/code/pointcloudfilter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/PointCloudFilter.dir/main.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PointCloudFilter.dir/main.cc.o -c /home/vboxuser/code/pointcloudfilter/main.cc

CMakeFiles/PointCloudFilter.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PointCloudFilter.dir/main.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vboxuser/code/pointcloudfilter/main.cc > CMakeFiles/PointCloudFilter.dir/main.cc.i

CMakeFiles/PointCloudFilter.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PointCloudFilter.dir/main.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vboxuser/code/pointcloudfilter/main.cc -o CMakeFiles/PointCloudFilter.dir/main.cc.s

# Object files for target PointCloudFilter
PointCloudFilter_OBJECTS = \
"CMakeFiles/PointCloudFilter.dir/main.cc.o"

# External object files for target PointCloudFilter
PointCloudFilter_EXTERNAL_OBJECTS =

PointCloudFilter: CMakeFiles/PointCloudFilter.dir/main.cc.o
PointCloudFilter: CMakeFiles/PointCloudFilter.dir/build.make
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_apps.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_outofcore.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_people.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libboost_system.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libboost_regex.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libqhull.so
PointCloudFilter: /usr/lib/libOpenNI.so
PointCloudFilter: /usr/lib/libOpenNI2.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkChartsCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkInfovisCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libfreetype.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libz.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libjpeg.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpng.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libtiff.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libexpat.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOGeometry-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOLegacy-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOPLY-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingLOD-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkViewsContext2D-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkViewsCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingContextOpenGL2-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingOpenGL2-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libflann_cpp.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_surface.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_keypoints.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_tracking.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_recognition.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_registration.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_stereo.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_segmentation.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_features.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_filters.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_sample_consensus.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_ml.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_visualization.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_search.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_kdtree.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_io.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_octree.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libpcl_common.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkInteractionWidgets-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersModeling-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkInteractionStyle-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersExtraction-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersStatistics-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkImagingFourier-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkalglib-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersHybrid-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkImagingGeneral-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkImagingSources-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkImagingHybrid-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingAnnotation-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkImagingColor-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingVolume-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOXML-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOXMLParser-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingContext2D-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeType-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libfreetype.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkImagingCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkRenderingCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonColor-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeometry-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersSources-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneral-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonComputationalGeometry-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkFiltersCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkIOImage-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonExecutionModel-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonDataModel-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonTransforms-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonMisc-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonMath-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonSystem-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkCommonCore-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtksys-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkDICOMParser-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libvtkmetaio-7.1.so.7.1p.1
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libz.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libGLEW.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libSM.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libICE.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libX11.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libXext.so
PointCloudFilter: /usr/lib/x86_64-linux-gnu/libXt.so
PointCloudFilter: CMakeFiles/PointCloudFilter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vboxuser/code/pointcloudfilter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable PointCloudFilter"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PointCloudFilter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PointCloudFilter.dir/build: PointCloudFilter

.PHONY : CMakeFiles/PointCloudFilter.dir/build

CMakeFiles/PointCloudFilter.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PointCloudFilter.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PointCloudFilter.dir/clean

CMakeFiles/PointCloudFilter.dir/depend:
	cd /home/vboxuser/code/pointcloudfilter/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vboxuser/code/pointcloudfilter /home/vboxuser/code/pointcloudfilter /home/vboxuser/code/pointcloudfilter/build /home/vboxuser/code/pointcloudfilter/build /home/vboxuser/code/pointcloudfilter/build/CMakeFiles/PointCloudFilter.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PointCloudFilter.dir/depend
