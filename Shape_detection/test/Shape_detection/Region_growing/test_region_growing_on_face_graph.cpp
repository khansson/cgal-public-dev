
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
bool test_region_growing_on_face_graph(int argc, char *argv[]) {

    using FT      = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;

    using Face_graph = CGAL::Surface_mesh<Point_3>;
    using Face_range = typename Face_graph::Face_range;

    using Connectivity   = CGAL::Shape_detection::Connectivity_on_face_graph<Face_graph>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Kernel, Face_graph>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Connectivity, Conditions>;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/cube.off");
    CGAL::set_ascii_mode(in);

    Face_graph face_graph;
    in >> face_graph;
    
    in.close();
    const Face_range face_range = face_graph.faces();

    CGAL_assertion(face_range.size() == 6);
    if (face_range.size() != 6) return false;

    // Create connectivity and conditions.
    using Vertex_to_point_map = typename Conditions::Vertex_to_point_map;
    const Vertex_to_point_map vertex_to_point_map(get(CGAL::vertex_point, face_graph));

    Connectivity connectivity(face_graph);
    Conditions     conditions(face_graph, max_distance_to_plane, normal_threshold, min_region_size, vertex_to_point_map);

    // Run region growing.
    Region_growing region_growing(face_range, connectivity, conditions);
    region_growing.detect();

    // Test data.
    const auto &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 6 && regions.size() == region_growing.number_of_regions());
    if (regions.size() != 6) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.are_valid(*region)) return false;

    const auto &unassigned_faces = region_growing.unassigned_items();

    CGAL_assertion(unassigned_faces.size() == 0 && unassigned_faces.size() == region_growing.number_of_unassigned_items());
    if (unassigned_faces.size() != 0) return false;

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
