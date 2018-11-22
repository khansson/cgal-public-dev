
// STL includes.
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

// Type declarations.
using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT                = Kernel::FT;
using Point_2           = Kernel::Point_2;
using Vector_2          = Kernel::Vector_2;
using Point_with_normal = std::pair<Point_2, Vector_2>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Pwn_vector        = std::vector<Point_with_normal>;
using Input_range       = CGAL::Iterator_range<typename Pwn_vector::iterator>;

using Traits            = CGAL::Shape_detection::Region_growing_traits<Input_range, Point_map, Kernel>;
using Connectivity      = CGAL::Shape_detection::Fuzzy_sphere_connectivity_on_points<Traits>;
using Conditions        = CGAL::Shape_detection::Propagation_conditions_on_points_2<Traits, Normal_map>;
using Region_growing    = CGAL::Shape_detection::Region_growing<Traits, Connectivity, Conditions>;
using Regions           = typename Region_growing::Region_range;

// Todo:
// What about including Iterator_range inside the class so that users do not need to wrap Pwn_vector by hand?
// Add ply exporter with colors.

int main(int argc, char *argv[]) {
    
    std::cout << std::endl << "region_growing_on_points_2 example started" << std::endl << std::endl;
    std::cout << "Note: if 0 points are loaded, please specify the path to the file data/points_2.pwn by hand!" << std::endl << std::endl;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/points_2.pwn");
    CGAL::set_ascii_mode(in);

    Pwn_vector pwn;
    double a,b,c,d,e,f;

    while (in >> a >> b >> c >> d >> e >> f)
        pwn.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
    
    Input_range input_range(pwn.begin(), pwn.end());
    std::cout << "* loaded " << pwn.size() << " points with normals" << std::endl;

    // Create instances of the classes Connectivity and Conditions
    const FT     search_radius        = FT(5);
    const FT     max_distance_to_line = FT(45) / FT(10);
    const FT     normal_thershold     = FT(7)  / FT(10);
    const size_t min_region_size      = 5;

    Connectivity connectivity(input_range, search_radius);
    Conditions   conditions(max_distance_to_line, normal_thershold, min_region_size);
    
    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.find_regions();

    // Print the number of regions found.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    std::cout << std::endl << "region_growing_on_points_2 example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
