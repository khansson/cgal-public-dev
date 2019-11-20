#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.

// General.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Points_2 = std::vector<Point_2>;

// Solver.
using MatrixFT = Eigen::SparseMatrix<FT>;
using Solver   = Eigen::SimplicialLDLT<MatrixFT>;

// Coordinates.
using Domain   = CGAL::Barycentric_coordinates::Delaunay_domain_2<
  Points_2, Kernel>;
using Harmonic = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<
  Points_2, Domain, Solver, Kernel>;

int main() {  

  // Construct a unit square.
  const Points_2 square = { 
    Point_2(0.0, 0.0), Point_2(1.0, 0.0),
    Point_2(1.0, 1.0), Point_2(0.0, 1.0) };

  // Instantiate a Delaunay domain.
  std::list<Point_2> list_of_seeds;
  list_of_seeds.push_back(Point_2(0.5, 0.5));
  
  Domain domain(square);
  domain.create(0.1, list_of_seeds);

  // Instantiate a sparse solver.
  Solver solver;

  // Compute harmonic coordinates at the vertices of the domain.
  Harmonic harmonic(
    square, domain, solver);
  harmonic.compute();
  
  // Create an std::vector to store coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(square.size());

  // Output harmonic coordinates.
  std::cout << std::endl << 
    "harmonic coordinates: " 
  << std::endl << std::endl;
  for (std::size_t k = 0; k < domain.number_of_vertices(); ++k) {
    coordinates.clear();
    harmonic.coordinates(
      k, std::back_inserter(coordinates));

    for (std::size_t i = 0; i < coordinates.size() - 1; ++i)
      std::cout << coordinates[i] << ", ";
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }

  // Evaluate harmonic coordinates at the barycenters of all finite elements
  // and output them one by one.
  std::cout << std::endl << 
    "harmonic coordinates evaluated at: " 
  << std::endl << std::endl;

  Points_2 barycenters;
  domain.get_barycenters(barycenters);
  for (const auto& barycenter : barycenters) {
    coordinates.clear();
    CGAL::Barycentric_coordinates::analytic_coordinates_2(
      square, barycenter, harmonic, std::back_inserter(coordinates));

    std::cout << barycenter << ": ";
    for (std::size_t i = 0; i < coordinates.size() - 1; ++i)
      std::cout << coordinates[i] << ", ";
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }
  std::cout << std::endl;
}
