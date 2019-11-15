#include <list>
#include <vector>
#include <string>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_weights_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

struct Info {
  
  Info(const std::string _name) : 
  name(_name) { }
  std::string name;
};

using Vertex_with_info = std::pair<Point_2, Info>;
using Vertices_2       = std::list<Vertex_with_info>;

using Point_with_id = std::pair<Point_2, std::size_t>;
using Points_2      = std::list<Point_with_id>;

using Vertex_map = CGAL::First_of_pair_property_map<Vertex_with_info>;
using Point_map  = CGAL::First_of_pair_property_map<Point_with_id>;

using Discrete_harmonic = 
  CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2<Vertices_2, Kernel, Vertex_map>;
using Policy = CGAL::Barycentric_coordinates::Computation_policy;

int main() {
  
  Kernel kernel;
  Point_map point_map;
  Vertex_map vertex_map;

  // Construct a unit square.
  const Vertices_2 square = {
    std::make_pair(Point_2(0, 0), Info("1")), std::make_pair(Point_2(1, 0), Info("2")), 
    std::make_pair(Point_2(1, 1), Info("3")), std::make_pair(Point_2(0, 1), Info("4"))
  };

  // Instantiate the class with discrete harmonic weights.
  // We do not check for edge cases since we know the exact positions 
  // of all our points.
  const Policy policy = Policy::PRECISE_COMPUTATION;
  Discrete_harmonic discrete_harmonic(square, policy, vertex_map);

  // Instantiate the center point of the unit square.
  const Point_2 center(FT(1) / FT(2), FT(1) / FT(2));

  // Compute discrete harmonic weights for the center point.
  std::list<FT> weights;
  discrete_harmonic(center, std::back_inserter(weights));

  std::cout << std::endl << "discrete harmonic weights (center): ";
  for (const FT weight : weights) std::cout << weight << " ";
  std::cout << std::endl;

  // Compute discrete harmonic coordinates for the center point.
  std::list<FT> coordinates;
  CGAL::Barycentric_coordinates::analytic_coordinates_2(
    square, center, discrete_harmonic, std::back_inserter(coordinates));

  std::cout << std::endl << "discrete harmonic coordinates (center): ";
  for (const FT coordinate : coordinates) std::cout << coordinate << " ";
  std::cout << std::endl;
  
  // Instantiate several interior points.
  const Points_2 interior_points = { 
    std::make_pair(Point_2(FT(1) / FT(5), FT(1) / FT(5)), 0), 
    std::make_pair(Point_2(FT(4) / FT(5), FT(1) / FT(5)), 1),
    std::make_pair(Point_2(FT(4) / FT(5), FT(4) / FT(5)), 2),
    std::make_pair(Point_2(FT(1) / FT(5), FT(4) / FT(5)), 3) };

  // Compute discrete harmonic weights for all the interior points at once.
  std::vector< std::vector<FT> > ws;
  ws.reserve(interior_points.size());
  CGAL::Barycentric_coordinates::analytic_weights_2(
    square, interior_points, discrete_harmonic, 
    std::back_inserter(ws), kernel, point_map);

  std::cout << std::endl << 
    "discrete harmonic weights (interior): " 
  << std::endl << std::endl;
  for (const auto& w : ws) {
    for (std::size_t i = 0; i < w.size() - 1; ++i) 
      std::cout << w[i] << ", ";
    std::cout << w[w.size() - 1] << std::endl;
  }

  // Compute discrete harmonic coordinates for all the interior points at once.
  std::vector< std::vector<FT> > bs;
  bs.reserve(interior_points.size());
  CGAL::Barycentric_coordinates::analytic_coordinates_2(
    square, interior_points, discrete_harmonic, 
    std::back_inserter(bs), kernel, point_map);

  std::cout << std::endl << 
    "discrete harmonic coordinates (interior): " 
  << std::endl << std::endl;
  for (const auto& b : bs) {
    for (std::size_t i = 0; i < b.size() - 1; ++i) 
      std::cout << b[i] << ", ";
    std::cout << b[b.size() - 1] << std::endl;
  }

  // Instantiate 2 boundary points on the second and fourth edges.
  const Point_2 e2(1, FT(4) / FT(5));
  const Point_2 e4(0, FT(4) / FT(5));

  // Compute discrete harmonic coordinates = boundary coordinates 
  // for these 2 points one by one.
  coordinates.clear();
  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    square, e2, std::back_inserter(coordinates), vertex_map);
  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    square, e4, std::back_inserter(coordinates), vertex_map);

  std::cout << std::endl << "boundary coordinates edge2 edge4: ";
  for (const FT coordinate : coordinates)
    std::cout << coordinate << " ";
  std::cout << std::endl;

  // Instantiate 6 other boundary points: 2 on the first and third edges respectively
  // and 4 at the vertices.
  const Points_2 es13 = {
    std::make_pair(Point_2(FT(1) / FT(2), 0), 1), // edges
    std::make_pair(Point_2(FT(1) / FT(2), 1), 3),

    std::make_pair(Point_2(0, 0), -1), // vertices
    std::make_pair(Point_2(1, 0), -1),
    std::make_pair(Point_2(1, 1), -1),
    std::make_pair(Point_2(0, 1), -1)
  };

  // Compute discrete harmonic coordinates = boundary coordinates 
  // for all these 6 points at once.
  bs.clear();
  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    square, es13, std::back_inserter(bs), kernel, vertex_map, point_map);

  std::cout << std::endl << 
    "boundary coordinates edge1 edge3 + vertices: " 
  << std::endl << std::endl;
  for (const auto& b : bs) {
    for (std::size_t i = 0; i < b.size() - 1; ++i) 
      std::cout << b[i] << ", ";
    std::cout << b[b.size() - 1] << std::endl;
  }

  // Instantiate 2 points outside the unit square - one from the left and one from the right.
  // Even if discrete harmonic coordinates may not be valid for some exterior points,
  // we can still do it.
  const Point_2 l(FT(-1) / FT(2), FT(1) / FT(2));
  const Point_2 r(FT(3)  / FT(2), FT(1) / FT(2));

  // Compute discrete harmonic coordinates for all the exterior points one by one.
  coordinates.clear();
  auto result = CGAL::Barycentric_coordinates::analytic_coordinates_2(
    square, l, discrete_harmonic, std::back_inserter(coordinates));
  result = CGAL::Barycentric_coordinates::analytic_coordinates_2(
    square, r, discrete_harmonic, std::back_inserter(coordinates));

  std::cout << std::endl << "discrete harmonic coordinates (exterior): ";
  for (const FT coordinate : coordinates) std::cout << coordinate << " ";
  std::cout << std::endl;

  // Return status of the last computation.
  const std::string status = (result ? "SUCCESS." : "FAILURE.");
  std::cout << std::endl << 
    "status of the last computation: " 
  << status << std::endl << std::endl;

  return EXIT_SUCCESS;
}
