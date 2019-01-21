// STL includes.
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Iterator_range.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_face_graph.h>

#include <CGAL/Shape_detection/Region_growing/Property_maps/Generic_index_to_item_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Property_maps/Hashable_item_to_index_property_map.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

#define USE_SURFACE_MESH

#if defined(USE_SURFACE_MESH)

    using Face_graph   = CGAL::Surface_mesh<Point_3>;
    using Face_range   = typename Face_graph::Face_range;

    using Connectivity = CGAL::Shape_detection::Connectivity_on_face_graph<Face_graph>;
    using Conditions   = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Kernel, Face_graph>;

#else

    using Face_graph        = CGAL::Polyhedron_3<Kernel>;
    using Face_range        = typename CGAL::Iterator_range<typename boost::graph_traits<Face_graph>::face_iterator>;
    using Index_to_face_map = CGAL::Shape_detection::Generic_index_to_item_property_map<Face_range>;
    using Face_to_index_map = CGAL::Shape_detection::Hashable_item_to_index_property_map<Face_range>;

    using Connectivity = CGAL::Shape_detection::Connectivity_on_face_graph<Face_graph, Face_range, Index_to_face_map, Face_to_index_map>;
    using Conditions   = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Kernel, Face_graph, Face_range, Index_to_face_map>;

#endif

using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Connectivity, Conditions>;

int main(int argc, char *argv[]) {

    std::cout << std::endl << "region_growing_on_face_graph example started" << std::endl << std::endl;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/face_graph.off");
    CGAL::set_ascii_mode(in);

    Face_graph face_graph;
    in >> face_graph;
    
    in.close();
    const Face_range face_range = CGAL::faces(face_graph);

    std::cout << "* cube mesh with " << face_range.size() << " faces is loaded" << std::endl;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1) / FT(10);
    const FT          normal_threshold      = FT(9) / FT(10);
    const std::size_t min_region_size       = 1;

    // Create instances of the classes Connectivity and Conditions.
    #if defined(USE_SURFACE_MESH)

        using Vertex_to_point_map = typename Conditions::Vertex_to_point_map;
        const Vertex_to_point_map vertex_to_point_map(get(CGAL::vertex_point, face_graph));

        Connectivity connectivity(face_graph);
        Conditions     conditions(face_graph, max_distance_to_plane, normal_threshold, min_region_size, vertex_to_point_map);

    #else

        const Index_to_face_map index_to_face_map(face_range);
        const Face_to_index_map face_to_index_map(face_range);
        
        Connectivity connectivity(face_graph, index_to_face_map, face_to_index_map);
        Conditions     conditions(face_graph, index_to_face_map, max_distance_to_plane, normal_threshold, min_region_size);

    #endif

    // Create an instance of the region growing class.
    Region_growing region_growing(face_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.detect();

    // Print the number of found regions.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    std::cout << std::endl << "region_growing_on_face_graph example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
