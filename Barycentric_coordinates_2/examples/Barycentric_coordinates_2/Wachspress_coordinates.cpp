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
  const std::size_t num_queries = 100;

  // Create vectors to store query points and polygon vertices.
  Points_2 queries, convex;

  // Generate a set of random query points.
  queries.reserve(num_queries);
  Generator generator(1.0);
  std::copy_n(
    generator, num_queries, std::back_inserter(queries));

  // Find the convex hull of the generated query points.
  // This convex hull gives the vertices of a convex polygon 
  // that contains all the generated points.
  CGAL::convex_hull_2(
    queries.begin(), queries.end(), std::back_inserter(convex));
  const std::size_t num_vertices = convex.size();

  // Instantiate the class with Wachspress weights.
  Wachspress wachspress(convex);
    
  // Compute Wachspress coordinates for all query points at once.
  std::vector< std::vector<FT> > bs;
  bs.reserve(queries.size());
  CGAL::Barycentric_coordinates::analytic_coordinates_2(
    convex, queries, wachspress, 
    std::back_inserter(bs), Kernel(), Point_map());

  // Output Wachspress coordinates.
  std::cout << std::endl << 
    "Wachspress coordinates (interior + boundary): " 
  << std::endl << std::endl;
  for (const auto& b : bs) {
    for (std::size_t i = 0; i < b.size() - 1; ++i)
      std::cout << b[i] << ", ";
    std::cout << b[b.size() - 1] << std::endl;
  }

  return EXIT_SUCCESS;
}
