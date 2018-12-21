
// STL includes.
#include <string>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Iterator_range.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_face_graph.h>

#include <CGAL/Shape_detection/Region_growing/Property_maps/Generic_index_to_item_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Property_maps/Hashable_item_to_index_property_map.h>

template<class Kernel>
bool test_region_growing_on_polyhedron(int argc, char *argv[]) {

    using FT      = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;

    using Polyhedron        = CGAL::Polyhedron_3<Kernel>;
    using Face_range        = typename CGAL::Iterator_range<typename boost::graph_traits<Polyhedron>::face_iterator>;
    using Index_to_face_map = CGAL::Shape_detection::Generic_index_to_item_property_map<Face_range>;
    using Face_to_index_map = CGAL::Shape_detection::Hashable_item_to_index_property_map<Face_range>;

    using Connectivity   = CGAL::Shape_detection::Connectivity_on_face_graph<Polyhedron, Face_range, Index_to_face_map, Face_to_index_map>;
    using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Kernel, Polyhedron, Face_range, Index_to_face_map>;
    using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Connectivity, Conditions>;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/polyhedron.off");
    CGAL::set_ascii_mode(in);

    Polyhedron polyhedron;
    in >> polyhedron;
    
    in.close();
    const Face_range face_range = CGAL::faces(polyhedron);

    CGAL_assertion(face_range.size() == 6);
    if (face_range.size() != 6) return false;

    // Create connectivity and conditions.
    const Index_to_face_map index_to_face_map(face_range);
    const Face_to_index_map face_to_index_map(face_range);
        
    Connectivity connectivity(polyhedron, index_to_face_map, face_to_index_map);
    Conditions     conditions(polyhedron, index_to_face_map, max_distance_to_plane, normal_threshold, min_region_size);

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
    if (!test_region_growing_on_polyhedron< CGAL::Simple_cartesian<double> >(argc, argv)) cartesian_double_test_success = false;
    
    std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
    CGAL_assertion(cartesian_double_test_success);

    // ------>

    bool exact_inexact_test_success = true;
    if (!test_region_growing_on_polyhedron<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)) exact_inexact_test_success = false;
    
    std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
    CGAL_assertion(exact_inexact_test_success);

    // ------>

    bool exact_exact_test_success = true;
    if (!test_region_growing_on_polyhedron<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv)) exact_exact_test_success = false;
    
    std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
    CGAL_assertion(exact_exact_test_success);
}
