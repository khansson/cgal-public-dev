
// STL includes.
#include <string>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_face_graph.h>

template<class Kernel>
bool test_region_growing_on_face_graph(int argc, char *argv[]) {

    using FT      = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;

    using Surface_mesh        = CGAL::Surface_mesh<Point_3>;
    using Input_range         = typename Surface_mesh::Face_range;
    using Face_descriptor     = typename Surface_mesh::Face_index;
    using Face_descriptor_map = CGAL::Identity_property_map<Face_descriptor>;

    using Connectivity   = CGAL::Shape_detection::Connectivity_on_face_graph<Surface_mesh, Face_descriptor_map>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Surface_mesh, Face_descriptor_map, Kernel>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Connectivity, Conditions>;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/cube.off");
    CGAL::set_ascii_mode(in);

    Surface_mesh surface_mesh;
    in >> surface_mesh;
    
    in.close();
    const Input_range input_range = surface_mesh.faces();

    CGAL_assertion(input_range.size() == 6);
    if (input_range.size() != 6) return false;

    // Create connectivity and conditions.
    Connectivity connectivity(surface_mesh);
    Conditions     conditions(surface_mesh, max_distance_to_plane, normal_threshold, min_region_size);

    // Run region growing.
    Region_growing region_growing(input_range, connectivity, conditions);
    region_growing.detect();

    // Test data.
    const auto &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 6 && regions.size() == region_growing.number_of_regions());
    if (regions.size() != 6) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.are_valid(*region)) return false;

    const auto &unclassified_faces = region_growing.unclassified_items();

    CGAL_assertion(unclassified_faces.size() == 0 && unclassified_faces.size() == region_growing.number_of_unclassified_items());
    if (unclassified_faces.size() != 0) return false;

    return true;
}

int main(int argc, char *argv[]) {

    // ------>

    bool cartesian_double_test_success = true;
    if (!test_region_growing_on_face_graph< CGAL::Simple_cartesian<double> >(argc, argv)) cartesian_double_test_success = false;
    
    std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
    CGAL_assertion(cartesian_double_test_success);

    // ------>

    bool exact_inexact_test_success = true;
    if (!test_region_growing_on_face_graph<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) exact_inexact_test_success = false;
    
    std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
    CGAL_assertion(exact_inexact_test_success);

    // ------>

    bool exact_exact_test_success = true;
    if (!test_region_growing_on_face_graph<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv)) exact_exact_test_success = false;
    
    std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
    CGAL_assertion(exact_exact_test_success);
}
