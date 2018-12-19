
// STL includes.
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_face_graph.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Face_range   = typename Surface_mesh::Face_range;
using Face_index   = typename Surface_mesh::Face_index;

using Connectivity   = CGAL::Shape_detection::Connectivity_on_face_graph<Surface_mesh>;
using Conditions     = CGAL::Shape_detection::Propagation_conditions_on_face_graph<Kernel, Surface_mesh>;
using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Connectivity, Conditions>;

using Color = CGAL::Color;

int main(int argc, char *argv[]) {

    std::cout << std::endl << "region_growing_on_surface_mesh example started" << std::endl << std::endl;

    // Load data.
    std::ifstream in(argc > 1 ? argv[1] : "../data/surface_mesh.off");
    CGAL::set_ascii_mode(in);

    Surface_mesh surface_mesh;
    in >> surface_mesh;
    
    in.close();
    const Face_range face_range = CGAL::faces(surface_mesh);

    std::cout << "* surface mesh with " << face_range.size() << " faces is loaded" << std::endl;

    // Default parameter values for the cube mesh below.
    const FT          max_distance_to_plane = FT(1);
    const FT          normal_threshold      = FT(7) / FT(10);
    const std::size_t min_region_size       = 5;

    // Create instances of the classes Connectivity and Conditions.
    using Vertex_to_point_map = typename Conditions::Vertex_to_point_map;
    const Vertex_to_point_map vertex_to_point_map(get(CGAL::vertex_point, surface_mesh));

    Connectivity connectivity(surface_mesh);
    Conditions     conditions(surface_mesh, max_distance_to_plane, normal_threshold, min_region_size, vertex_to_point_map);

    // Create an instance of the region growing class.
    Region_growing region_growing(face_range, connectivity, conditions);

    // Run the algorithm.
    region_growing.detect();

    // Print the number of found regions.
    std::cerr << "* " << region_growing.number_of_regions() << " regions have been found" << std::endl;

    // Get all found regions.
    const auto &regions = region_growing.regions();
    
    // Save result to a file in the user-provided path if any.
    srand(time(NULL));
    if (argc > 2) {
        
        bool created;
        typename Surface_mesh::template Property_map<Face_index, Color> face_color;
        boost::tie(face_color, created) = surface_mesh.template add_property_map<Face_index, Color>("f:color", Color(0, 0, 0));

        if (!created) {
            
            std::cout << std::endl << "region_growing_on_surface_mesh example finished" << std::endl << std::endl;
            EXIT_FAILURE;
        }

        const std::string path     = argv[2];
        const std::string fullpath = path + "regions.off";

        std::ofstream out(fullpath);

        // Iterate through all regions.
        for (auto region = regions.begin(); region != regions.end(); ++region) {
            
            // Generate a random color.
            const Color color(
                static_cast<unsigned char>(rand() % 256),
                static_cast<unsigned char>(rand() % 256),
                static_cast<unsigned char>(rand() % 256));

            // Iterate through all region items.
            for (auto index : *region)
                face_color[Face_index(index)] = color;
        }

        out << surface_mesh;
        out.close();

        std::cerr << "* surface mesh is saved in " << fullpath << std::endl;
    }

    std::cout << std::endl << "region_growing_on_surface_mesh example finished" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
