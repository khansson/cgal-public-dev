
// STL includes.
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_surface_mesh.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using Mesh         = CGAL::Surface_mesh<Point_3>;
using Vertex_index = typename Mesh::Vertex_index;
using Face_index   = typename Mesh::Face_index;
using Face_range   = typename Mesh::Face_range;
using Faces        = std::vector<Face_index>;
using Input_range  = CGAL::Iterator_range<typename Faces::iterator>;
using Element_map  = CGAL::Identity_property_map<Face_index>;

using Traits         = CGAL::Shape_detection::Region_growing_traits<Input_range, Element_map, Kernel>;
using Connectivity   = CGAL::Shape_detection::Connectivity_on_surface_mesh<Traits, Mesh>;
using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_surface_mesh<Traits, Mesh>;
using Region_growing = CGAL::Shape_detection::Region_growing<Traits, Connectivity, Conditions>;

// Todo:
// What about including Iterator_range inside the class so that users do not need to wrap Pwn_vector by hand?

int main(int argc, char *argv[]) {

    std::cout << std::endl << "region_growing_on_surface_mesh example started" << std::endl << std::endl;

    // Create simple mesh of a cube.
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

    std::cout << "* cube mesh with " << mesh.faces().size() << " faces and " << mesh.vertices().size() << " vertices is created" << std::endl;

    // Set up a list of faces as the input range.
    Faces cube_faces;
    const Face_range &mesh_faces = mesh.faces();
    
    for (auto it = mesh_faces.begin(); it != mesh_faces.end(); ++it) cube_faces.push_back(*it);
    Input_range input_range(cube_faces.begin(), cube_faces.end());

    std::cout << "* input range of faces is set" << std::endl;

    // Create instances of the classes Connectivity and Conditions.
    const FT max_distance_to_plane = FT(1) / FT(10);
    const FT normal_threshold      = FT(9) / FT(10);

    Connectivity connectivity(mesh);
    Conditions   conditions(mesh, max_distance_to_plane, normal_threshold);

    // Create an instance of the region growing class.
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.find_regions();

    // Print the number of regions found.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    std::cout << std::endl << "region_growing_on_surface_mesh example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
