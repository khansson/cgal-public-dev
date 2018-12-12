
// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

/*
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h> */

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

// Todo:
/* // should be fixed!
using Input_range = CGAL::Point_set_3<Point_3>;
using Point_map   = typename Input_range::Point_map;
using Normal_map  = typename Input_range::Vector_map; */

using Point_with_normal = std::pair<Point_3, Vector_3>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Traits         = CGAL::Shape_detection::Region_growing_traits<Input_range, Point_map, Kernel>;
using Connectivity   = CGAL::Shape_detection::Nearest_neighbor_connectivity_on_points<Traits>;
using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_points_3<Traits, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Traits, Connectivity, Conditions>;
using Regions        = typename Region_growing::Region_range;

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
    
        Output_rep (const ::Color &c) : c(c) { }

        std::ostream& operator() (std::ostream &out) const {
            
            if (is_ascii(out)) out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
            else out.write(reinterpret_cast<const char*>(&c), sizeof(c));
            
            return out;
        }
    };

} // namespace CGAL

int main(int argc, char *argv[]) {
    
    std::cout << std::endl << "region_growing_on_points_3 example started" << std::endl << std::endl;
    std::cout << "Note: if 0 points are loaded, please specify the path to the file data/points_3.xyz by hand!" << std::endl << std::endl;

    // Load xyz data either from a local folder or a user-provided file.
    std::ifstream in(argc > 1 ? argv[1] : "../data/points_3.xyz");
    CGAL::set_ascii_mode(in);

    /*
    const bool with_normal_map = true;
    Input_range input_range(with_normal_map);

    in >> input_range;    
    in.close(); */

    Input_range input_range;
    Point_3 point; Vector_3 normal;

    while (in >> point >> normal)
        input_range.push_back(std::make_pair(point, normal));
    
    in.close();
    std::cout << "* loaded " << input_range.size() << " points with normals" << std::endl;

    // Create instances of the classes Connectivity and Conditions.
    const size_t num_neighbors         = 100;
    const FT     max_distance_to_plane = FT(5) / FT(10);
    const FT     normal_threshold      = FT(9) / FT(10);
    const size_t min_region_size       = 3;

    Connectivity connectivity(input_range, num_neighbors);
    Conditions   conditions(max_distance_to_plane, normal_threshold, min_region_size);

    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.find_regions();

    // Print the number of found regions.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    // Get the list of found regions.
    const Regions &regions = region_growing.regions();

    Pwc_vector pwc;
    srand(time(NULL));

    // Iterate through all regions.
    for (auto region = regions.begin(); region != regions.end(); ++region) {
        
        // Generate a random color.
        const Color color = CGAL::make_array(
                                static_cast<unsigned char>(rand() % 256),
                                static_cast<unsigned char>(rand() % 256),
                                static_cast<unsigned char>(rand() % 256));

        // Iterate through all region elements.
        for (auto element : *region) {
            
            const Point_3 &point = get(Point_map(), element);
            pwc.push_back(std::make_pair(point, color));
        }
    }

    // Save result to a file in the user-provided path if any.
    if (argc > 2) {
        
        const std::string path     = argv[2];
        const std::string fullpath = path + "regions_3.ply";

        std::ofstream out(fullpath);

        CGAL::set_ascii_mode(out);
        CGAL::write_ply_points_with_properties(out, pwc,
                                        CGAL::make_ply_point_writer(PLY_Point_map()),
                                        std::make_tuple(PLY_Color_map(), 
                                        CGAL::PLY_property<unsigned char>("red"),
                                        CGAL::PLY_property<unsigned char>("green"),
                                        CGAL::PLY_property<unsigned char>("blue")));

        std::cerr << "* found regions are saved in " << fullpath << std::endl;
    }

    std::cout << std::endl << "region_growing_on_points_3 example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
