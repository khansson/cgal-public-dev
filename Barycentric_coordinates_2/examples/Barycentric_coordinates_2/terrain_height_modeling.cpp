#include <map>
#include <vector>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_weights_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.

// General.
using Kernel     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Projection = CGAL::Projection_traits_xy_3<Kernel>;

using FT    = typename Projection::FT;
using Point = typename Projection::Point_2;

// Triangulation.
using FB  = CGAL::Delaunay_mesh_face_base_2<Projection>;
using VB  = CGAL::Triangulation_vertex_base_2<Projection>;
using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<Projection, TDS>;

using Vertex_handle = typename CDT::Vertex_handle;

// Mesher.
using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
using Mesher   = CGAL::Delaunay_mesher_2<CDT, Criteria>;

// Coordinates.
using Points     = std::vector<Point>;
using Mean_value = CGAL::Barycentric_coordinates::Mean_value_weights_2<Points, Projection>;

// Interpolation.
using Interpolation_traits   = CGAL::Interpolation_traits_2<Projection>;
using Vertex_function_value  = std::map<Point, FT, typename Projection::Less_xy_2>;
using Function_value_access  = CGAL::Data_access<Vertex_function_value>;
using Point_with_coordinate  = std::pair<Point, FT>;
using Points_with_coordinate = std::vector<Point_with_coordinate>;

int main() {

  // Construct a polygon that bounds a three-dimensional terrain.
  // Note that z-coordinate of each vertex represents the height function.
  // Projection in 2D is performed automatically by the Projection traits class.
  const Points polygon = {
    Point(0.03, 0.05, 0.000), Point(0.07, 0.04, 10.00), Point(0.10, 0.04, 20.00),
    Point(0.14, 0.04, 30.00), Point(0.17, 0.07, 40.00), Point(0.19, 0.09, 50.00),
    Point(0.22, 0.11, 60.00), Point(0.25, 0.11, 70.00), Point(0.27, 0.10, 80.00),
    Point(0.30, 0.07, 90.00), Point(0.31, 0.04, 100.0), Point(0.34, 0.03, 110.0),
    Point(0.37, 0.02, 120.0), Point(0.40, 0.03, 130.0), Point(0.42, 0.04, 140.0),
    Point(0.44, 0.07, 150.0), Point(0.45, 0.10, 160.0), Point(0.46, 0.13, 170.0),
    Point(0.46, 0.19, 180.0), Point(0.47, 0.26, 190.0), Point(0.47, 0.31, 200.0),
    Point(0.47, 0.35, 210.0), Point(0.45, 0.37, 220.0), Point(0.41, 0.38, 230.0),
    Point(0.38, 0.37, 240.0), Point(0.35, 0.36, 250.0), Point(0.32, 0.35, 260.0),
    Point(0.30, 0.37, 270.0), Point(0.28, 0.39, 280.0), Point(0.25, 0.40, 290.0),
    Point(0.23, 0.39, 300.0), Point(0.21, 0.37, 310.0), Point(0.21, 0.34, 320.0),
    Point(0.23, 0.32, 330.0), Point(0.24, 0.29, 340.0), Point(0.27, 0.24, 350.0),
    Point(0.29, 0.21, 360.0), Point(0.29, 0.18, 370.0), Point(0.26, 0.16, 380.0),
    Point(0.24, 0.17, 390.0), Point(0.23, 0.19, 400.0), Point(0.24, 0.22, 410.0),
    Point(0.24, 0.25, 420.0), Point(0.21, 0.26, 430.0), Point(0.17, 0.26, 440.0),
    Point(0.12, 0.24, 450.0), Point(0.07, 0.20, 460.0), Point(0.03, 0.15, 470.0),
    Point(0.01, 0.10, 480.0), Point(0.02, 0.07, 490.0)
  };

  // Create a Delaunay triangulation with polygon edges as constraints.
  CDT cdt;
  std::vector<Vertex_handle> vhs;
  vhs.reserve(polygon.size());
  for (const auto& vertex : polygon) 
    vhs.push_back(cdt.insert(vertex));

  for(std::size_t i = 0; i < vhs.size(); ++i) {
    const std::size_t ip = (i + 1) % vhs.size();
    cdt.insert_constraint(vhs[i], vhs[ip]);
  }

  // Refine this triangulation.
  Mesher mesher(cdt);
  mesher.set_criteria(Criteria(0.01, 0.01));
  mesher.refine_mesh();

  // Associate each polygon vertex with the corresponding function value.
  Vertex_function_value vertex_function_value;
  for(const auto& vertex : polygon)
    vertex_function_value.insert(
      std::make_pair(vertex, vertex.z()));

  Points_with_coordinate boundary;
  boundary.resize(polygon.size());

  // Store all generated interior points with interpolated data here.
  Points queries;
  queries.reserve(cdt.number_of_vertices());

  // Instantiate the class with mean value weights.
  Mean_value mean_value(polygon);

  // Compute mean value coordinates and use them to interpolate data 
  // from the polygon's boundary to its interior.
  std::vector<FT> coordinates;
  coordinates.reserve(polygon.size());
  for (auto vh = cdt.finite_vertices_begin();
  vh != cdt.finite_vertices_end(); ++vh) {
    
    coordinates.clear();
    const auto& query = vh->point();
    CGAL::Barycentric_coordinates::analytic_coordinates_2(
    polygon, query, mean_value, std::back_inserter(coordinates), Projection());

    for (std::size_t i = 0; i < polygon.size(); ++i) 
      boundary[i] = std::make_pair(polygon[i], coordinates[i]);
    
    const FT f = CGAL::linear_interpolation(
      boundary.begin(), boundary.end(), FT(1), 
      Function_value_access(vertex_function_value));
    
    queries.push_back(Point(query.x(), query.y(), f));
  }
  
  // Output interpolated heights.
  std::cout << std::endl << "interpolated heights: " << std::endl << std::endl;
  for (const auto& query : queries)
    std::cout << query.z() << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
