// STL includes.
#include <string>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Surface_mesh.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_polygon_mesh(int argc, char *argv[]) {

  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  using Surface_mesh = CGAL::Surface_mesh<Point_3>;
  using Face_range   = typename Surface_mesh::Face_range;

  using Connectivity   = SD::Polygon_mesh_adjacent_faces_connectivity<Surface_mesh>;
  using Conditions     = SD::Polygon_mesh_least_squares_plane_fit_conditions<Kernel, Surface_mesh>;
  using Region_growing = SD::Region_growing<Face_range, Connectivity, Conditions>;

  // Default parameter values for the data file polygon_mesh.off.
  const FT          max_distance_to_plane = FT(1);
  const FT          normal_threshold      = FT(7) / FT(10);
  const std::size_t min_region_size       = 5;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "../data/polygon_mesh.off");
  CGAL::set_ascii_mode(in);

  Surface_mesh surface_mesh;
  in >> surface_mesh;
    
  in.close();
  const Face_range face_range = CGAL::faces(surface_mesh);

  CGAL_assertion(face_range.size() == 32245);
  if (face_range.size() != 32245) 
    return false;

  // Create connectivity and conditions.
  Connectivity connectivity(
    surface_mesh);

  using Vertex_to_point_map = typename Conditions::Vertex_to_point_map;
  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, surface_mesh));

  Conditions conditions(
    surface_mesh, 
    max_distance_to_plane, normal_threshold, min_region_size, 
    vertex_to_point_map);

  // Run region growing.
  Region_growing region_growing(face_range, connectivity, conditions);
  region_growing.detect();

  // Test data.
  const auto& regions = region_growing.regions();

  CGAL_assertion(
    regions.size() == 334 && 
    regions.size() == region_growing.number_of_regions());

  if (regions.size() != 334) 
    return false;

  for (auto region = regions.begin(); region != regions.end(); ++region)
    if (!conditions.is_valid_region(*region)) 
      return false;

  const auto& unassigned_faces = region_growing.unassigned_items();

  CGAL_assertion(
    unassigned_faces.size() == 858 && 
    unassigned_faces.size() == region_growing.number_of_unassigned_items());

  if (unassigned_faces.size() != 858) 
    return false;

  return true;
}

int main(int argc, char *argv[]) {

  // ------>

  bool cartesian_double_test_success = true;
  if (!test_region_growing_on_polygon_mesh< CGAL::Simple_cartesian<double> >(argc, argv)) 
    cartesian_double_test_success = false;
    
  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
  CGAL_assertion(cartesian_double_test_success);

  // ------>

  bool exact_inexact_test_success = true;
  if (!test_region_growing_on_polygon_mesh<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) 
    exact_inexact_test_success = false;
    
  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
  CGAL_assertion(exact_inexact_test_success);
}
