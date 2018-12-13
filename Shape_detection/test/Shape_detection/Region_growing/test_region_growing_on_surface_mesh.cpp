
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
#include <CGAL/Shape_detection/Region_growing/Region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_surface_mesh.h>

template<class Kernel>
bool test_region_growing_on_surface_mesh() {

    using FT      = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;

    using Mesh           = CGAL::Surface_mesh<Point_3>;
    using Input_range    = typename Mesh::Face_range;
    using Vertex_index   = typename Mesh::Vertex_index;
    using Face_index     = typename Mesh::Face_index;
    using Face_index_map = CGAL::Identity_property_map<Face_index>;

    using Traits         = CGAL::Shape_detection::Region_growing_traits<Input_range, Face_index_map, Kernel>;
    using Connectivity   = CGAL::Shape_detection::Connectivity_on_surface_mesh<Traits, Mesh>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_surface_mesh<Traits, Mesh>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Traits, Connectivity, Conditions>;
    using Regions        = typename Region_growing::Region_range;

    // Default parameter values for the cube mesh below.
    const FT max_distance_to_plane = FT(1) / FT(10);
    const FT normal_threshold      = FT(9) / FT(10);
    
    // Create a cube mesh.
    Mesh mesh;

    const Vertex_index v0 = mesh.add_vertex(Point_3(0, 0, 0));
    const Vertex_index v1 = mesh.add_vertex(Point_3(0, 1, 0));
    const Vertex_index v2 = mesh.add_vertex(Point_3(1, 0, 0));
    const Vertex_index v3 = mesh.add_vertex(Point_3(1, 1, 0));
    const Vertex_index v4 = mesh.add_vertex(Point_3(0, 0, 1));
    const Vertex_index v5 = mesh.add_vertex(Point_3(0, 1, 1));
    const Vertex_index v6 = mesh.add_vertex(Point_3(1, 0, 1));
    const Vertex_index v7 = mesh.add_vertex(Point_3(1, 1, 1));

    mesh.add_face(v0, v3, v2);
    mesh.add_face(v0, v3, v1);
    mesh.add_face(v0, v6, v2);
    mesh.add_face(v0, v6, v4);
    mesh.add_face(v0, v5, v1);
    mesh.add_face(v0, v5, v4);
    mesh.add_face(v7, v2, v3);
    mesh.add_face(v7, v2, v6);
    mesh.add_face(v7, v1, v3);
    mesh.add_face(v7, v1, v5);
    mesh.add_face(v7, v4, v5);
    mesh.add_face(v7, v4, v6);

    const Input_range &input_range = mesh.faces();

    CGAL_assertion(input_range.size() == 6);
    if (input_range.size() != 6) return false;

    // Create connectivity and conditions.
    Connectivity connectivity(mesh);
    Conditions   conditions(mesh, max_distance_to_plane, normal_threshold);

    // Run region growing.
    Region_growing region_growing(input_range, connectivity, conditions);
    region_growing.detect();

    const Regions &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 6);
    if (regions.size() != 6) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.are_valid(*region)) return false;
    return true;
}

int main(int argc, char *argv[]) {

    // ------>

    bool cartesian_double_test_success = true;
    if (!test_region_growing_on_surface_mesh< CGAL::Simple_cartesian<double> >()) cartesian_double_test_success = false;
    
    std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
    CGAL_assertion(cartesian_double_test_success);

    // ------>

    bool exact_inexact_test_success = true;
    if (!test_region_growing_on_surface_mesh<CGAL::Exact_predicates_inexact_constructions_kernel>()) exact_inexact_test_success = false;
    
    std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
    CGAL_assertion(exact_inexact_test_success);

    // ------>

    bool exact_exact_test_success = true;
    if (!test_region_growing_on_surface_mesh<CGAL::Exact_predicates_exact_constructions_kernel>()) exact_exact_test_success = false;
    
    std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
    CGAL_assertion(exact_exact_test_success);
}
