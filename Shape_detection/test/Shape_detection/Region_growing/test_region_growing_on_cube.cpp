// STL includes.
#include <string>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_vector.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_cube(int argc, char *argv[]) {

    using FT      = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;

    using Polyhedron        = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>;
    using Face_range        = typename CGAL::Iterator_range<typename boost::graph_traits<Polyhedron>::face_iterator>;

    using Connectivity   = SD::Polygon_mesh_adjacent_faces_connectivity<Polyhedron, Face_range>;
    using Conditions     = SD::Polygon_mesh_least_squares_plane_fit_conditions<Kernel, Polyhedron, Face_range>;
    using Region_growing = SD::Region_growing<Face_range, Connectivity, Conditions>;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/cube.off");
    CGAL::set_ascii_mode(in);

    Polyhedron polyhedron;
    in >> polyhedron;
    
    in.close();
    const Face_range face_range = CGAL::faces(polyhedron);

    CGAL_assertion(face_range.size() == 6);
    if (face_range.size() != 6) return false;

    // Create connectivity and conditions.
    Connectivity connectivity(polyhedron);
    Conditions     conditions(polyhedron, max_distance_to_plane, normal_threshold, min_region_size);

    // Run region growing.
    Region_growing region_growing(face_range, connectivity, conditions);
    region_growing.detect();

    // Test data.
    const auto &regions = region_growing.regions();

    CGAL_assertion(regions.size() == 6 && regions.size() == region_growing.number_of_regions());
    if (regions.size() != 6) return false;

    for (auto region = regions.begin(); region != regions.end(); ++region)
        if (!conditions.is_valid_region(*region)) return false;

    const auto &unassigned_faces = region_growing.unassigned_items();

    CGAL_assertion(unassigned_faces.size() == 0 && unassigned_faces.size() == region_growing.number_of_unassigned_items());
    if (unassigned_faces.size() != 0) return false;

    return true;
}

int main(int argc, char *argv[]) {

    // ------>

    bool cartesian_double_test_success = true;
    if (!test_region_growing_on_cube< CGAL::Simple_cartesian<double> >(argc, argv)) cartesian_double_test_success = false;
    
    std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
    CGAL_assertion(cartesian_double_test_success);

    // ------>

    bool exact_inexact_test_success = true;
    if (!test_region_growing_on_cube<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) exact_inexact_test_success = false;
    
    std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
    CGAL_assertion(exact_inexact_test_success);

    // ------>

    bool exact_exact_test_success = true;
    if (!test_region_growing_on_cube<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv)) exact_exact_test_success = false;
    
    std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
    CGAL_assertion(exact_exact_test_success);
}
