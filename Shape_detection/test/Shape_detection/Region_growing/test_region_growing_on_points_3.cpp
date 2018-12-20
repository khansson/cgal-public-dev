
// STL includes.
#include <string>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

template<class Kernel>
bool test_region_growing_on_points_3(int argc, char *argv[]) {

    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;

    using Input_range = CGAL::Point_set_3<Point_3>;
    using Point_map   = typename Input_range::Point_map;
    using Normal_map  = typename Input_range::Vector_map;

    using Connectivity   = CGAL::Shape_detection::Nearest_neighbor_connectivity_on_points<Kernel, Input_range, Point_map>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_points_3<Kernel, Input_range, Point_map, Normal_map>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Connectivity, Conditions>;

    // Default parameter values for the data file points_3.xyz.
    const std::size_t num_neighbors         = 100;
    const FT          max_distance_to_plane = FT(5) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 3;
    
    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/points_3.xyz");
    CGAL::set_ascii_mode(in);

    const bool with_normal_map = true;
    Input_range input_range(with_normal_map);

    in >> input_range;
    in.close();

    CGAL_assertion(input_range.size() == 300000);
    if (input_range.size() != 300000) return false;

    // Create connectivity and conditions.
    Connectivity connectivity(input_range, num_neighbors, input_range.point_map());
    Conditions     conditions(input_range, max_distance_to_plane, normal_threshold, min_region_size, input_range.point_map(), input_range.normal_map());

    // Run region growing.
    Region_growing region_growing(input_range, connectivity, conditions);
    region_growing.detect();

    // Test data.
    const auto &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 1108);
    if (regions.size() != 1108) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.are_valid(*region)) return false;

    const auto &unassigned_points = region_growing.unassigned_items();
 
    CGAL_assertion(unassigned_points.size() == 1063 && unassigned_points.size() == region_growing.number_of_unassigned_items());
    if (unassigned_points.size() != 1063) return false;

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
}
