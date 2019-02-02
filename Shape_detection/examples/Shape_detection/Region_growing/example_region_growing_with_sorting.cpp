// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Point_3  = typename Kernel::Point_3;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Connectivity   = SD::Points_fuzzy_sphere_connectivity<Kernel, Input_range, Point_map>;
using Conditions     = SD::Points_2_least_squares_line_fit_conditions<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = SD::Points_2_least_squares_line_fit_sorting<Kernel, Input_range, Point_map, Connectivity>;
using Region_growing = SD::Region_growing<Input_range, Connectivity, Conditions, typename Sorting::Seed_map>;

using Color            = CGAL::cpp11::array<unsigned char, 3>;
using Point_with_color = std::pair<Point_3, Color>;
using Pwc_vector       = std::vector<Point_with_color>;
using PLY_Point_map    = CGAL::First_of_pair_property_map<Point_with_color>;
using PLY_Color_map    = CGAL::Second_of_pair_property_map<Point_with_color>;

// Define how a color should be stored.
namespace CGAL {
    
  template<class F>
  struct Output_rep< ::Color, F > {
        
    const ::Color& c;
    static const bool is_specialized = true;
    
    Output_rep(const ::Color& c) : c(c) { }

    std::ostream& operator()(std::ostream& out) const {
            
      if (is_ascii(out)) 
        out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
      else 
        out.write(reinterpret_cast<const char*>(&c), sizeof(c));
            
      return out;
    }
  };

} // namespace CGAL

int main(int argc, char *argv[]) {
    
  std::cout << std::endl << 
    "region_growing_with_sorting example started" 
  << std::endl << std::endl;
    
  std::cout << 
    "Note: if 0 points are loaded, please specify the path to the file data/points_2.xyz by hand!" 
  << std::endl << std::endl;

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "../data/points_2.xyz");
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

  Sorting sorting(input_range, connectivity);
  sorting.sort();
    
  // Create an instance of the region growing class.
  Region_growing region_growing(input_range, connectivity, conditions, sorting.seed_map());

  // Run the algorithm.
  region_growing.detect();

  // Print the number of found regions.
  std::cerr << "* " << region_growing.number_of_regions() << 
    " regions have been found" 
  << std::endl;

  // Get all found regions.
  const auto& regions = region_growing.regions();

  Pwc_vector pwc;
  srand(time(NULL));

  // Iterate through all regions.
  for (auto region = regions.begin(); region != regions.end(); ++region) {
        
    // Generate a random color.
    const Color color = 
      CGAL::make_array(
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256));

    // Iterate through all region items.
    for (auto index : *region) {
      const auto& key = *(input_range.begin() + index);

      const Point_2& point = get(Point_map(), key);
      pwc.push_back(std::make_pair(Point_3(point.x(), point.y(), 0), color));
    }
  }

  // Save result to a file in the user-provided path if any.
  if (argc > 2) {
        
    const std::string path     = argv[2];
    const std::string fullpath = path + "regions_points_2.ply";

    std::ofstream out(fullpath);

    CGAL::set_ascii_mode(out);
    CGAL::write_ply_points_with_properties(
      out, pwc,
      CGAL::make_ply_point_writer(PLY_Point_map()),
        std::make_tuple(
          PLY_Color_map(), 
          CGAL::PLY_property<unsigned char>("red"),
          CGAL::PLY_property<unsigned char>("green"),
          CGAL::PLY_property<unsigned char>("blue")));

    std::cout << 
      "* found regions are saved in " 
    << fullpath << std::endl;
    out.close();
  }

  std::cout << std::endl << 
    "region_growing_on_points_2 example finished" 
  << std::endl << std::endl;
  
  return EXIT_SUCCESS;
}
