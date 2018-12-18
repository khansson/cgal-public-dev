
// STL includes.
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_face_graph.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Generic_index_to_item_property_map.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

// #define POLYHEDRON_GRAPH_TYPE

#if defined(POLYHEDRON_GRAPH_TYPE)
    
    using Face_graph      = CGAL::Polyhedron_3<Kernel>;
    using Input_range     = typename CGAL::Iterator_range< boost::graph_traits<Face_graph>::face_iterator>;
    using Face_descriptor = typename Face_graph::Face_handle;

#else

    using Face_graph      = CGAL::Surface_mesh<Point_3>;
    using Input_range     = typename Face_graph::Face_range;
    using Face_descriptor = typename Face_graph::Face_index;

#endif

using Face_descriptor_map = CGAL::Identity_property_map<Face_descriptor>;
using Index_to_item_map   = CGAL::Shape_detection::Generic_index_to_item_property_map<Input_range>;

using Connectivity   = CGAL::Shape_detection::Connectivity_on_face_graph<Face_graph, Face_descriptor_map, Input_range, Index_to_item_map>;
using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Face_graph, Face_descriptor_map, Kernel, Input_range, Index_to_item_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Connectivity, Conditions>;

int main(int argc, char *argv[]) {

    std::cout << std::endl << "region_growing_on_face_graph example started" << std::endl << std::endl;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/cube.off");
    CGAL::set_ascii_mode(in);

    Face_graph face_graph;
    in >> face_graph;
    
    in.close();
    const Input_range input_range = CGAL::faces(face_graph);

    std::cout << "* cube mesh with " << input_range.size() << " faces is loaded" << std::endl;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Create instances of the classes Connectivity and Conditions.
    const Index_to_item_map index_to_item_map(input_range);

    Connectivity connectivity(face_graph, index_to_item_map);
    Conditions     conditions(face_graph, index_to_item_map, max_distance_to_plane, normal_threshold, min_region_size);

    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.detect();

    // Print the number of found regions.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    std::cout << std::endl << "region_growing_on_face_graph example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
