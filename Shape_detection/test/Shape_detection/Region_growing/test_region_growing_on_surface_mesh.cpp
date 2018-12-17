
// STL includes.
#include <string>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Surface_mesh.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_face_graph.h>

template<class Kernel>
bool test_region_growing_on_surface_mesh(int argc, char *argv[]) {

    using FT      = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;

    using Mesh           = CGAL::Surface_mesh<Point_3>;
    using Input_range    = typename Mesh::Face_range;
    using Vertex_index   = typename Mesh::Vertex_index;
    using Face_index     = typename Mesh::Face_index;
    using Face_index_map = CGAL::Identity_property_map<Face_index>;

    using Connectivity   = CGAL::Shape_detection::Connectivity_on_face_graph<Mesh, Face_index_map>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Mesh, Face_index_map, Kernel>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Connectivity, Conditions>;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/cube.off");
    CGAL::set_ascii_mode(in);

    Mesh mesh;
    in >> mesh;
    
    in.close();
    const Input_range &input_range = mesh.faces();

    CGAL_assertion(input_range.size() == 6);
    if (input_range.size() != 6) return false;

    // Create connectivity and conditions.
    Connectivity connectivity(mesh);
    Conditions   conditions(mesh, max_distance_to_plane, normal_threshold, min_region_size);

    // Run region growing.
    Region_growing region_growing(input_range, connectivity, conditions);
    region_growing.detect();

    const auto &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 6);
    if (regions.size() != 6) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.are_valid(*region)) return false;
    return true;
}

int main(int argc, char *argv[]) {

    // ------>

    bool cartesian_double_test_success = true;
    if (!test_region_growing_on_surface_mesh< CGAL::Simple_cartesian<double> >(argc, argv)) cartesian_double_test_success = false;
    
    std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
    CGAL_assertion(cartesian_double_test_success);

    // ------>

    bool exact_inexact_test_success = true;
    if (!test_region_growing_on_surface_mesh<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) exact_inexact_test_success = false;
    
    std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
    CGAL_assertion(exact_inexact_test_success);

    // ------>

    bool exact_exact_test_success = true;
    if (!test_region_growing_on_surface_mesh<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv)) exact_exact_test_success = false;
    
    std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
    CGAL_assertion(exact_exact_test_success);
}
