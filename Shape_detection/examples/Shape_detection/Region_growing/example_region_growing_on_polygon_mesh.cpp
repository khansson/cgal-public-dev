// STL includes.
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/memory.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Color = CGAL::Color;

// Choose the type of a container for a polygon mesh.
#define USE_SURFACE_MESH

#if defined(USE_SURFACE_MESH)

    using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
    using Face_range   = typename Polygon_mesh::Face_range;

    using Connectivity = SD::Polygon_mesh_adjacent_faces_connectivity<Polygon_mesh>;
    using Conditions   = SD::Polygon_mesh_least_squares_plane_fit_conditions<Kernel, Polygon_mesh>;

#else

    using Polygon_mesh = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>;
    using Face_range   = typename CGAL::Iterator_range<typename boost::graph_traits<Polygon_mesh>::face_iterator>;
    
    using Connectivity = SD::Polygon_mesh_adjacent_faces_connectivity<Polygon_mesh, Face_range>;
    using Conditions   = SD::Polygon_mesh_least_squares_plane_fit_conditions<Kernel, Polygon_mesh, Face_range>;

#endif

using Region_growing = SD::Region_growing<Face_range, Connectivity, Conditions>;

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "region_growing_on_polygon_mesh example started" 
  << std::endl << std::endl;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "../data/polygon_mesh.off");
  CGAL::set_ascii_mode(in);

  Polygon_mesh polygon_mesh;
  in >> polygon_mesh;
    
  in.close();
  const Face_range face_range = CGAL::faces(polygon_mesh);

  std::cout << 
    "* polygon mesh with " 
  << face_range.size() << 
    " faces is loaded" 
  << std::endl;

  // Default parameter values for the data file polygon_mesh.off.
  const FT          max_distance_to_plane = FT(1);
  const FT          normal_threshold      = FT(7) / FT(10);
  const std::size_t min_region_size       = 5;

  // Create instances of the classes Connectivity and Conditions.
  Connectivity connectivity(
    polygon_mesh);

  using Vertex_to_point_map = typename Conditions::Vertex_to_point_map;
  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, polygon_mesh));

  Conditions conditions(
    polygon_mesh, 
    max_distance_to_plane, normal_threshold, min_region_size, 
    vertex_to_point_map);

  // Create an instance of the region growing class.
  Region_growing region_growing(face_range, connectivity, conditions);

  // Run the algorithm.
  region_growing.detect();

  // Print the number of found regions.
  std::cout << "* " << region_growing.number_of_regions() << 
    " regions have been found" 
  << std::endl;

  // Save result in a file only if it is stored in CGAL::Surface_mesh.
  #if defined(USE_SURFACE_MESH)
  
    using Face_index = typename Polygon_mesh::Face_index;

    // Get all found regions.
    const auto& regions = region_growing.regions();
      
    // Save result to a file in the user-provided path if any.
    srand(time(NULL));
    if (argc > 2) {

      bool created;
      typename Polygon_mesh::template Property_map<Face_index, Color> face_color;
      boost::tie(face_color, created) = 
        polygon_mesh.template add_property_map<Face_index, Color>(
          "f:color", Color(0, 0, 0));

      if (!created) {
              
        std::cout << std::endl << 
          "region_growing_on_polygon_mesh example finished" 
        << std::endl << std::endl;
          
        return EXIT_FAILURE;
      }

      const std::string path     = argv[2];
      const std::string fullpath = path + "regions_polygon_mesh.off";

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

      out << polygon_mesh;
      out.close();

      std::cout << 
        "* polygon mesh is saved in " 
      << fullpath << std::endl;
    }

  #endif

  std::cout << std::endl << 
    "region_growing_on_polygon_mesh example finished" 
  << std::endl << std::endl;
    
  return EXIT_SUCCESS;
}
