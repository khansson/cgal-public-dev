# Created by the script cgal_create_cmake_script This is the CMake script for
# compiling a CGAL application.

project(Efficient_RANSAC_Tests)

cmake_minimum_required(VERSION 3.1)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}"
     VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "test_efficient_RANSAC_cone_connected_component.cpp")
  create_single_source_cgal_program(
    "test_efficient_RANSAC_cone_parameters.cpp")

  create_single_source_cgal_program(
    "test_efficient_RANSAC_cylinder_connected_component.cpp")
  create_single_source_cgal_program(
    "test_efficient_RANSAC_cylinder_parameters.cpp")

  create_single_source_cgal_program(
    "test_efficient_RANSAC_plane_connected_component.cpp")
  create_single_source_cgal_program(
    "test_efficient_RANSAC_plane_parameters.cpp")

  create_single_source_cgal_program(
    "test_efficient_RANSAC_sphere_connected_component.cpp")
  create_single_source_cgal_program(
    "test_efficient_RANSAC_sphere_parameters.cpp")

  create_single_source_cgal_program(
    "test_efficient_RANSAC_torus_connected_component.cpp")
  create_single_source_cgal_program(
    "test_efficient_RANSAC_torus_parameters.cpp")

  create_single_source_cgal_program(
    "test_efficient_RANSAC_scene.cpp")
  create_single_source_cgal_program(
    "test_efficient_RANSAC_regularization.cpp")

else()
  message(
    STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
