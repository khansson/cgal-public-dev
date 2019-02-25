// STL includes.
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;

using Input_range = CGAL::Point_set_3<Point_3>;
using Point_map   = typename Input_range::Point_map;
using Normal_map  = typename Input_range::Vector_map;

using Connectivity   = SD::Points_k_nearest_neighbors_connectivity<Kernel, Input_range, Point_map>;
using Conditions     = SD::Points_3_least_squares_plane_fit_conditions<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = SD::Points_3_least_squares_plane_fit_sorting<Kernel, Input_range, Connectivity, Point_map>;
using Region_growing = SD::Region_growing<Input_range, Connectivity, Conditions, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "../data/points_3.xyz");
  CGAL::set_ascii_mode(in);

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);

  in >> input_range;
  in.close();

  // Default parameter values for the data file points_3.xyz.
  const size_t num_neighbors         = 100;
  const FT     max_distance_to_plane = FT(5) / FT(10);
  const FT     normal_threshold      = FT(9) / FT(10);
  const size_t min_region_size       = 3;

  // Create instances of the classes Connectivity and Conditions.
  Connectivity connectivity(
    input_range, 
    num_neighbors, 
    input_range.point_map());
    
  Conditions conditions(
    input_range, 
    max_distance_to_plane, normal_threshold, min_region_size, 
    input_range.point_map(), input_range.normal_map());

  // Sort point indices.
  Sorting sorting(
    input_range, connectivity,
    input_range.point_map());
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, connectivity, conditions, 
    sorting.seed_map());

  // Run the algorithm.
  std::vector<typename Region_growing::Item_indices> regions;
  region_growing.detect(std::back_inserter(regions));

  region_growing.release_memory();
  CGAL_assertion(regions.size() == 1180);

  const bool exact_inexact_test_success = (regions.size() == 1180);
  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;

  return EXIT_SUCCESS;
}
