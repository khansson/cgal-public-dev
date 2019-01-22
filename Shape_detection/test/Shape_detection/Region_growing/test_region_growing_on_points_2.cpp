// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_points_2(int argc, char *argv[]) {
    
  using FT       = typename Kernel::FT;
  using Point_2  = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;

  using Point_with_normal = std::pair<Point_2, Vector_2>;
  using Input_range       = std::vector<Point_with_normal>;
  using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
  using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

  using Connectivity   = SD::Points_fuzzy_sphere_connectivity<Kernel, Input_range, Point_map>;
  using Conditions     = SD::Points_2_least_squares_line_fit_conditions<Kernel, Input_range, Point_map, Normal_map>;
  using Region_growing = SD::Region_growing<Input_range, Connectivity, Conditions>;

  // Default parameter values for the data file points_2.xyz.
  const FT          search_radius        = FT(5);
  const FT          max_distance_to_line = FT(45) / FT(10);
  const FT          normal_threshold     = FT(7)  / FT(10);
  const std::size_t min_region_size      = 5;
    
  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "../data/points_2.xyz");
  CGAL::set_ascii_mode(in);

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
    
  in.close();

  CGAL_assertion(input_range.size() == 3634);
  if (input_range.size() != 3634) 
    return false;

  // Create connectivity and conditions.
  Connectivity connectivity(
    input_range, 
    search_radius);
  
  Conditions conditions(
    input_range, 
    max_distance_to_line, normal_threshold, min_region_size);

  // Run region growing.
  Region_growing region_growing(input_range, connectivity, conditions);
  region_growing.detect();

  // Test data.
  const auto& regions = region_growing.regions();

  CGAL_assertion(
    regions.size() == 65 && 
    regions.size() == region_growing.number_of_regions());

  if (regions.size() != 65) 
    return false;

  for (auto region = regions.begin(); region != regions.end(); ++region)
    if (!conditions.is_valid_region(*region)) 
      return false;

  const auto& unassigned_points = region_growing.unassigned_items();

  CGAL_assertion(
    unassigned_points.size() == 82 && 
    unassigned_points.size() == region_growing.number_of_unassigned_items());

  if (unassigned_points.size() != 82) 
    return false;

  return true;
}

int main(int argc, char *argv[]) {

  // ------>

  bool cartesian_double_test_success = true;
  if (!test_region_growing_on_points_2< CGAL::Simple_cartesian<double> >(argc, argv)) 
    cartesian_double_test_success = false;
    
  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
  CGAL_assertion(cartesian_double_test_success);

  // ------>

  bool exact_inexact_test_success = true;
  if (!test_region_growing_on_points_2<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) 
    exact_inexact_test_success = false;
    
  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
  CGAL_assertion(exact_inexact_test_success);

  // ------>

  bool exact_exact_test_success = true;
  if (!test_region_growing_on_points_2<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv)) 
    exact_exact_test_success = false;
    
  std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
  CGAL_assertion(exact_exact_test_success);
}
