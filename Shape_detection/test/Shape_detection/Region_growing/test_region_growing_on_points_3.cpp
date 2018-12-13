
// STL includes.
#include <string>
#include <utility>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

template<class Kernel>
bool test_region_growing_on_points_3(int argc, char *argv[]) {

    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;
    using Vector_3 = typename Kernel::Vector_3;

    using Point_with_normal = std::pair<Point_3, Vector_3>;
    using Input_range       = std::vector<Point_with_normal>;
    using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
    using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

    using Traits         = CGAL::Shape_detection::Region_growing_traits<Input_range, Point_map, Kernel>;
    using Connectivity   = CGAL::Shape_detection::Nearest_neighbor_connectivity_on_points<Traits>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_points_3<Traits, Normal_map>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Traits, Connectivity, Conditions>;
    using Regions        = typename Region_growing::Region_range;

    // Default parameter values for the data file points_3.xyz.
    const size_t num_neighbors         = 100;
    const FT     max_distance_to_plane = FT(5) / FT(10);
    const FT     normal_threshold      = FT(9) / FT(10);
    const size_t min_region_size       = 3;
    
    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/points_3.xyz");
    CGAL::set_ascii_mode(in);

    Point_3  point;
    Vector_3 normal;

    Input_range input_range;
    while (in >> point >> normal) 
        input_range.push_back(std::make_pair(point, normal));
    
    in.close();

    CGAL_assertion(input_range.size() == 300000);
    if (input_range.size() != 300000) return false;

    // Create connectivity and conditions.
    Connectivity connectivity(input_range, num_neighbors);
    Conditions   conditions(max_distance_to_plane, normal_threshold, min_region_size);

    // Run region growing.
    Region_growing region_growing(input_range, connectivity, conditions);
    region_growing.detect();

    const Regions &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 1108);
    if (regions.size() != 1108) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.are_valid(*region)) return false;
    return true;
}

int main(int argc, char *argv[]) {

    // ------>

    bool cartesian_double_test_success = true;
    if (!test_region_growing_on_points_3< CGAL::Simple_cartesian<double> >(argc, argv)) cartesian_double_test_success = false;
    
    std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
    CGAL_assertion(cartesian_double_test_success);

    // ------>

    bool exact_inexact_test_success = true;
    if (!test_region_growing_on_points_3<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) exact_inexact_test_success = false;
    
    std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
    CGAL_assertion(exact_inexact_test_success);

    // ------>

    /* // super slow test!
    bool exact_exact_test_success = true;
    if (!test_region_growing_on_points_3<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv)) exact_exact_test_success = false;
    
    std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
    CGAL_assertion(exact_exact_test_success); */
}
