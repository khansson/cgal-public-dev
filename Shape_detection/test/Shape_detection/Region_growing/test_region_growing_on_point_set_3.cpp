// STL includes.
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/property_map.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_point_set_3(int argc, char *argv[]) {

  using FT       = typename Kernel::FT;
  using Point_3  = typename Kernel::Point_3;

  using Input_range = CGAL::Point_set_3<Point_3>;
  using Point_map   = typename Input_range::Point_map;
  using Normal_map  = typename Input_range::Vector_map;

  using Neighbor_query = SD::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
  using Region_type    = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
  using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type>;

  // Default parameter values for the data file point_set_3.xyz.
  const std::size_t k                     = 100;
  const FT          max_distance_to_plane = FT(5) / FT(10);
  const FT          angle_threshold       = FT(25);
  const std::size_t min_region_size       = 3;
    
  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "../data/point_set_3.xyz");
  CGAL::set_ascii_mode(in);

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);

  in >> input_range;
  in.close();

  CGAL_assertion(input_range.size() == 300000);
  if (input_range.size() != 300000) 
    return false;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range, 
    k, 
    input_range.point_map());

  Region_type region_type(
    input_range, 
    max_distance_to_plane, angle_threshold, min_region_size, 
    input_range.point_map(), input_range.normal_map());

  // Run region growing.
  Region_growing region_growing(input_range, neighbor_query, region_type);
  
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  // Test data.
  CGAL_assertion(regions.size() == 1190);
  if (regions.size() != 1190) 
    return false;

  for (const auto& region : regions)
    if (!region_type.is_valid_region(region)) 
      return false;

  std::vector<std::size_t> unassigned_points;
  region_growing.output_unassigned_items(std::back_inserter(unassigned_points));

  CGAL_assertion(unassigned_points.size() == 1134);
  if (unassigned_points.size() != 1134) 
    return false;

  return true;
}

int main(int argc, char *argv[]) {

  // ------>

  bool cartesian_double_test_success = true;
  if (!test_region_growing_on_point_set_3< CGAL::Simple_cartesian<double> >(argc, argv)) 
    cartesian_double_test_success = false;
    
  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
  CGAL_assertion(cartesian_double_test_success);

  // ------>

  bool exact_inexact_test_success = true;
  if (!test_region_growing_on_point_set_3<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) 
    exact_inexact_test_success = false;
    
  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
  CGAL_assertion(exact_inexact_test_success);
}
