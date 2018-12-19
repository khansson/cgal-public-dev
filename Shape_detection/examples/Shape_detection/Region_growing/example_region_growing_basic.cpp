
// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing.h>

// Type declarations.
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Connectivity   = CGAL::Shape_detection::Fuzzy_sphere_connectivity_on_points<Kernel, Input_range, Point_map>;
using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_points_2<Kernel, Input_range, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Connectivity, Conditions>;

int main(int argc, char *argv[]) {
    
    std::cout << std::endl << "region_growing_basic example started" << std::endl << std::endl;
    std::cout << "Note: if 0 points are loaded, please specify the path to the file data/points_2.xyz by hand!" << std::endl << std::endl;

    // Load xyz data either from a local folder or a user-provided file.
    std::ifstream in(argc > 1 ? argv[1] : "../data/points_2.xyz");
    CGAL::set_ascii_mode(in);

    Input_range input_range;
    FT a, b, c, d, e, f;

    while (in >> a >> b >> c >> d >> e >> f)
        input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
    
    in.close();
    std::cout << "* loaded " << input_range.size() << " points with normals" << std::endl;

    // Create instances of the classes Connectivity and Conditions.
    const FT     search_radius        = FT(5);
    const FT     max_distance_to_line = FT(45) / FT(10);
    const FT     normal_threshold     = FT(7)  / FT(10);
    const size_t min_region_size      = 5;

    Connectivity connectivity(input_range, search_radius);
    Conditions     conditions(input_range, max_distance_to_line, normal_threshold, min_region_size);
    
    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.detect();

    // Print the number of found regions.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    std::cout << std::endl << "region_growing_basic example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
