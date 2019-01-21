// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_points.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

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

using Timer = CGAL::Timer;

void benchmark_region_growing_on_points_2(const size_t test_count, const Input_range &input_range, 
const FT search_radius, const FT max_distance_to_line, const FT normal_threshold, const size_t min_region_size) {

    // Create instances of the classes Connectivity and Conditions.
    Connectivity connectivity(input_range, search_radius);
    Conditions   conditions(input_range, max_distance_to_line, normal_threshold, min_region_size);
    
    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    Timer timer;
    timer.start();
    region_growing.detect();
    timer.stop();

    // Compute the number of points assigned to all found regions.
    size_t number_of_assigned_points = 0;
    const auto &regions = region_growing.regions();

    for (auto region = regions.begin(); region != regions.end(); ++region)
        number_of_assigned_points += region->size();

    // Print statistics.
    std::cout << "Test #"                          << test_count                                  << std::endl;
    std::cout << "  search_radius = "              << search_radius                               << std::endl;
    std::cout << "  min_region_size = "            << min_region_size                             << std::endl;
    std::cout << "  max_distance_to_line = "       << max_distance_to_line                        << std::endl;
    std::cout << "  normal_threshold = "           << normal_threshold                            << std::endl;
    std::cout << "  -----"                                                                        << std::endl;
    std::cout << "  Time elapsed: "                << timer.time()                                << std::endl;
    std::cout << "  Number of detected regions: "  << region_growing.number_of_regions()          << std::endl;
    std::cout << "  Number of assigned points: "   << number_of_assigned_points                   << std::endl;
    std::cout << "  Number of unassigned points: " << region_growing.number_of_unassigned_items() << std::endl;
    std::cout << std::endl << std::endl;
}

int main(int argc, char *argv[]) {
    
    // Load xyz data either from a local folder or a user-provided file.
    std::ifstream in(argc > 1 ? argv[1] : "../data/points_2.xyz");
    CGAL::set_ascii_mode(in);

    Input_range input_range;
    FT a, b, c, d, e, f;

    while (in >> a >> b >> c >> d >> e >> f) 
        input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));

    in.close();

    // Default parameter values for the data file points_2.xyz.
    const FT     max_distance_to_line = FT(45) / FT(10);
    const FT     normal_threshold     = FT(7)  / FT(10);
    const size_t min_region_size      = 5;

    // Run benchmarks.
    benchmark_region_growing_on_points_2(1, input_range, FT(1), max_distance_to_line, normal_threshold, min_region_size);
    benchmark_region_growing_on_points_2(2, input_range, FT(3), max_distance_to_line, normal_threshold, min_region_size);
    benchmark_region_growing_on_points_2(3, input_range, FT(6), max_distance_to_line, normal_threshold, min_region_size);
    benchmark_region_growing_on_points_2(4, input_range, FT(9), max_distance_to_line, normal_threshold, min_region_size);
}
