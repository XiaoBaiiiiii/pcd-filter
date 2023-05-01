cmake_minimum_required(VERSION 3.10)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

project(PointCloudFilter)

set(CMAKE_CXX_STANDARD 14)

find_package(PCL REQUIRED)

include_directories(${PCL_INCLUDE_DIRS}, CMAKE_CURRENT_SOURCE_DIR)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(PointCloudFilter main.cc)
target_link_libraries(PointCloudFilter ${PCL_LIBRARIES})