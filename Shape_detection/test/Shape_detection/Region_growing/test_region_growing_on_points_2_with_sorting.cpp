// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Connectivity   = SD::Points_fuzzy_sphere_connectivity<Kernel, Input_range, Point_map>;
using Conditions     = SD::Points_2_least_squares_line_fit_conditions<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = SD::Points_2_least_squares_line_fit_sorting<Kernel, Input_range, Connectivity, Point_map>;
using Region_growing = SD::Region_growing<Input_range, Connectivity, Conditions, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "../data/points_2.xyz");
  CGAL::set_ascii_mode(in);

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
  in.close();

  // Default parameter values for the data file points_2.xyz.
  const FT     search_radius        = FT(5);
  const FT     max_distance_to_line = FT(45) / FT(10);
  const FT     normal_threshold     = FT(7)  / FT(10);
  const size_t min_region_size      = 5;

  // Create instances of the classes Connectivity and Conditions.
  Connectivity connectivity(
    input_range, 
    search_radius);

  Conditions conditions(
    input_range, 
    max_distance_to_line, normal_threshold, min_region_size);

  // Sort point indices.
  Sorting sorting(
    input_range, connectivity);
  sorting.sort();
    
  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, connectivity, conditions, 
    sorting.seed_map());

  // Run the algorithm.
  std::vector<typename Region_growing::Item_indices> regions;
  region_growing.detect(std::back_inserter(regions));

  region_growing.release_memory();
  CGAL_assertion(regions.size() == 64);

  const bool cartesian_double_test_success = (regions.size() == 64);
  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;

  return EXIT_SUCCESS;
}
