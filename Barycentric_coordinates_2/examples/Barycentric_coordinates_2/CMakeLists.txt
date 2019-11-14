# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Barycentric_coordinates_2_Examples)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "segment_coordinates.cpp")
  create_single_source_cgal_program(
    "triangle_coordinates.cpp")
  create_single_source_cgal_program(
      "wachspress_coordinates.cpp")

else()
  message(WARNING 
  "This program requires the CGAL library, and will not be compiled.")
endif()
