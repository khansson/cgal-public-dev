
// STL includes.
#include <string>
#include <cstdlib>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_surface_mesh.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

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

int main(int argc, char *argv[]) {

    std::cout << std::endl << "region_growing_on_surface_mesh example started" << std::endl << std::endl;

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
    std::cout << "* cube mesh with " << mesh.faces().size() << " faces and " << mesh.vertices().size() << " vertices is created" << std::endl;

    // Create instances of the classes Connectivity and Conditions.
    const FT max_distance_to_plane = FT(1) / FT(10);
    const FT normal_threshold      = FT(9) / FT(10);

    Connectivity connectivity(mesh);
    Conditions   conditions(mesh, max_distance_to_plane, normal_threshold);

    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.find_regions();

    // Print the number of found regions.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    std::cout << std::endl << "region_growing_on_surface_mesh example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
