#include <vector>
#include <CGAL/property_map.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_weights_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Points_2 = std::vector<Point_2>;

using Point_map  = CGAL::Identity_property_map<Point_2>;
using Creator    = CGAL::Creator_uniform_2<FT, Point_2>;
using Generator  = CGAL::Random_points_in_square_2<Point_2, Creator>;
using Wachspress = CGAL::Barycentric_coordinates::Wachspress_weights_2<Points_2, Kernel>;

int main() {
  
  // Choose how many random query points we want to generate.
  const std::size_t num_query_points = 1000;

  // Create vectors to store query points and polygon vertices.
  Points_2 query_points, polygon;

  // Generate a set of random query points.
  query_points.reserve(num_query_points);
  Generator generator(1.0);
  std::copy_n(
    generator, num_query_points, std::back_inserter(query_points));

  // Find the convex hull of the generated query points.
  // This convex hull gives the vertices of a convex polygon 
  // that contains all the generated points.
  CGAL::convex_hull_2(
    query_points.begin(), query_points.end(), std::back_inserter(polygon));
  const std::size_t num_vertices = polygon.size();

  // Instantiate the class with Wachspress weights.
  Wachspress wachspress(polygon);
    
  // Compute Wachspress coordinates for all queries at once.
  std::vector< std::vector<FT> > coordinates;
  coordinates.reserve(query_points.size());
  CGAL::Barycentric_coordinates::analytic_coordinates_2(
    polygon, query_points, wachspress, 
    std::back_inserter(coordinates), Kernel(), Point_map());

  // Output Wachspress coordinates.
  std::cout << std::endl << "Wachspress coordinates: " << std::endl << std::endl;
  for (const auto& vec : coordinates) {
    for (const FT coordinate : vec)
      std::cout << coordinate << ", ";
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
