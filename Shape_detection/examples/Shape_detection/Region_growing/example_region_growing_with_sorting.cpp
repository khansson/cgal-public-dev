// STL includes.
#include <list>
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

// Type declarations.
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Neighbor_query = CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_line_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = CGAL::Shape_detection::Point_set::Least_squares_line_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {
    
  std::cout << std::endl << 
    "region_growing_with_sorting example started" 
  << std::endl << std::endl;
    
  std::cout << 
    "Note: if 0 points are loaded, please specify the path to the file data/point_set_2.xyz by hand!" 
  << std::endl << std::endl;

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "../data/point_set_2.xyz");
  CGAL::set_ascii_mode(in);

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
    
  in.close();
  std::cout << 
    "* loaded " 
  << input_range.size() << 
    " points with normals" 
  << std::endl;

  // Default parameter values for the data file point_set_2.xyz.
  const FT     sphere_radius        = FT(5);
  const FT     max_distance_to_line = FT(45) / FT(10);
  const FT     max_accepted_angle   = FT(45);
  const size_t min_region_size      = 5;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    input_range, 
    sphere_radius);

  Region_type region_type(
    input_range, 
    max_distance_to_line, max_accepted_angle, min_region_size);

  // Sort point indices.
  Sorting sorting(
    input_range, neighbor_query);

  sorting.sort();
    
  // Create an instance of the region growing class
  // that is using the sorted seeding order.
  Region_growing region_growing(
    input_range, neighbor_query, region_type,
    sorting.seed_map());

  // Run the algorithm and store its results in 
  // the user-defined container outside the class.
  std::list< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  // Release all internal memory.
  region_growing.release_memory();

  // Print the number of found regions.
  std::cerr << "* " << regions.size() << 
    " regions have been found" 
  << std::endl;

  std::cout << std::endl << 
    "region_growing_with_sorting example finished" 
  << std::endl << std::endl;
  
  return EXIT_SUCCESS;
}
