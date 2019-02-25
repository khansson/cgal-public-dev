// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection::Point_set;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Neighbor_query = SD::Sphere_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = SD::Least_squares_line_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type>;

using Timer = CGAL::Timer;

void benchmark_region_growing_on_points_2(
  const size_t test_count, 
  const Input_range& input_range, 
  const FT search_radius, 
  const FT max_distance_to_line, 
  const FT normal_threshold, 
  const size_t min_region_size) {

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    input_range, 
    search_radius);

  Region_type region_type(
    input_range, 
    max_distance_to_line, 
    normal_threshold, 
    min_region_size);
    
  // Create an instance of the region growing class.
  Region_growing region_growing(input_range, neighbor_query, region_type);

  // Run the algorithm.
  Timer timer;
  timer.start();
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  timer.stop();

  // Compute the number of points assigned to all found regions.
  size_t number_of_assigned_points = 0;
  for (const auto& region : regions)
    number_of_assigned_points += region.size();

  std::vector<std::size_t> unassigned_items;
  region_growing.output_unassigned_items(std::back_inserter(unassigned_items));

  // Print statistics.
  std::cout << "Test #"                          << test_count                                  << std::endl;
  std::cout << "  search_radius = "              << search_radius                               << std::endl;
  std::cout << "  min_region_size = "            << min_region_size                             << std::endl;
  std::cout << "  max_distance_to_line = "       << max_distance_to_line                        << std::endl;
  std::cout << "  normal_threshold = "           << normal_threshold                            << std::endl;
  std::cout << "  -----"                                                                        << std::endl;
  std::cout << "  Time elapsed: "                << timer.time()                                << std::endl;
  std::cout << "  Number of detected regions: "  << regions.size()                              << std::endl;
  std::cout << "  Number of assigned points: "   << number_of_assigned_points                   << std::endl;
  std::cout << "  Number of unassigned points: " << unassigned_items.size()                     << std::endl;
  std::cout << std::endl << std::endl;
}

int main(int argc, char *argv[]) {
    
  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "../data/point_set_2.xyz");
  CGAL::set_ascii_mode(in);

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f) 
    input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));

  in.close();

  // Default parameter values for the data file point_set_2.xyz.
  const FT     max_distance_to_line = FT(45) / FT(10);
  const FT     normal_threshold     = FT(7)  / FT(10);
  const size_t min_region_size      = 5;

  // Run benchmarks.
  benchmark_region_growing_on_points_2(1, input_range, FT(1), 
  max_distance_to_line, normal_threshold, min_region_size);

  benchmark_region_growing_on_points_2(2, input_range, FT(3), 
  max_distance_to_line, normal_threshold, min_region_size);

  benchmark_region_growing_on_points_2(3, input_range, FT(6), 
  max_distance_to_line, normal_threshold, min_region_size);

  benchmark_region_growing_on_points_2(4, input_range, FT(9), 
  max_distance_to_line, normal_threshold, min_region_size);
}
