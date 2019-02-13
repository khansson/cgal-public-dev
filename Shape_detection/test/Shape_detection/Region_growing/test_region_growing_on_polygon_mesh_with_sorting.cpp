// STL includes.
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
using Face_range   = typename Polygon_mesh::Face_range;

using Connectivity   = SD::Polygon_mesh_adjacent_faces_connectivity<Polygon_mesh>;
using Conditions     = SD::Polygon_mesh_least_squares_plane_fit_conditions<Kernel, Polygon_mesh>;
using Sorting        = SD::Polygon_mesh_least_squares_plane_fit_sorting<Kernel, Polygon_mesh, Connectivity>;
using Region_growing = SD::Region_growing<Face_range, Connectivity, Conditions, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "../data/polygon_mesh.off");
  CGAL::set_ascii_mode(in);

  Polygon_mesh polygon_mesh;
  in >> polygon_mesh;
    
  in.close();
  const Face_range face_range = CGAL::faces(polygon_mesh);

  // Default parameter values for the data file polygon_mesh.off.
  const FT          max_distance_to_plane = FT(1);
  const FT          normal_threshold      = FT(7) / FT(10);
  const std::size_t min_region_size       = 5;

  // Create instances of the classes Connectivity and Conditions.
  Connectivity connectivity(polygon_mesh);

  using Vertex_to_point_map = typename Conditions::Vertex_to_point_map;
  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, polygon_mesh));

  Conditions conditions(
    polygon_mesh, 
    max_distance_to_plane, normal_threshold, min_region_size, 
    vertex_to_point_map);

  // Sort face indices.
  Sorting sorting(
    polygon_mesh, connectivity,
    vertex_to_point_map);
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    face_range, connectivity, conditions, 
    sorting.seed_map());

  // Run the algorithm.
  std::vector<typename Region_growing::Item_indices> regions;
  region_growing.detect(std::back_inserter(regions));

  region_growing.release_memory();
  CGAL_assertion(regions.size() == 321);

  const bool exact_exact_test_success = (regions.size() == 321);
  std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;

  return EXIT_SUCCESS;
}
