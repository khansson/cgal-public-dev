// STL includes.
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using Input_range = CGAL::Point_set_3<Point_3>;
using Point_map   = typename Input_range::Point_map;
using Normal_map  = typename Input_range::Vector_map;

using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query, Region_type>;

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
    "region_growing_on_point_set_3 example started" 
  << std::endl << std::endl;
    
  std::cout << 
    "Note: if 0 points are loaded, please specify the path to the file data/point_set_3.xyz by hand!" 
  << std::endl << std::endl;

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "../data/point_set_3.xyz");
  CGAL::set_ascii_mode(in);

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);

  in >> input_range;
  in.close();

  std::cout << 
    "* loaded " 
  << input_range.size() << 
    " points with normals" 
  << std::endl;

  // Default parameter values for the data file point_set_3.xyz.
  const size_t k                     = 100;
  const FT     max_distance_to_plane = FT(5) / FT(10);
  const FT     max_accepted_angle    = FT(25);
  const size_t min_region_size       = 3;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    input_range, 
    k, 
    input_range.point_map());
    
  Region_type region_type(
    input_range, 
    max_distance_to_plane, max_accepted_angle, min_region_size,
    input_range.point_map(), input_range.normal_map());

  // Create an instance of the region growing class.
  Region_growing region_growing(input_range, neighbor_query, region_type);

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  // Print the number of found regions.
  std::cout << "* " << regions.size() << 
    " regions have been found" 
  << std::endl;

  Pwc_vector pwc;
  srand(time(NULL));

  // Iterate through all regions.
  for (const auto& region : regions) {
        
    // Generate a random color.
    const Color color = 
      CGAL::make_array(
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256));

    // Iterate through all region items.
    for (const auto index : region) {
      const auto& key = *(input_range.begin() + index);
            
      const Point_3& point = get(input_range.point_map(), key);
      pwc.push_back(std::make_pair(point, color));
    }
  }

  // Save result to a file in the user-provided path if any.
  if (argc > 2) {
        
    const std::string path     = argv[2];
    const std::string fullpath = path + "regions_point_set_3.ply";

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

    std::cout << "* found regions are saved in " << fullpath << std::endl;
    out.close();
  }

  // Get all unassigned points.
  std::vector<std::size_t> unassigned_items;
  region_growing.output_unassigned_items(
    std::back_inserter(unassigned_items));

  // Print the number of unassigned points.
  std::cerr << "* " << unassigned_items.size() << 
    " points do not belong to any region" 
  << std::endl;

  // Store all unassigned points.
  std::vector<Point_3> unassigned_points;

  for (const auto index : unassigned_items) {
    const auto& key = *(input_range.begin() + index);

    const Point_3& point = get(input_range.point_map(), key);
    unassigned_points.push_back(point);
  }

  std::cout << "* " << unassigned_points.size() << 
    " unassigned points are stored" 
  << std::endl;
  
  std::cout << std::endl << 
    "region_growing_on_point_set_3 example finished" 
  << std::endl << std::endl;

  return EXIT_SUCCESS;
}
