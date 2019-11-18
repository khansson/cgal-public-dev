#include <vector>
#include <iterator>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Points_2        = std::vector<Point_2>;
using Point_map       = CGAL::Identity_property_map<Point_2>;
using Output_iterator = std::back_insert_iterator< std::vector<FT> >;

using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

int main() {

  Kernel traits;
  Point_map point_map;

  // Create a set of vertices.
  const Points_2 vertices = { 
    Point_2(0.0, 0.0), Point_2(0.75, 0.25), Point_2(0.5, 0.5), Point_2(0.4, -0.2) };

  // Create a set of query points.
  const Points_2 queries = { 
    Point_2(0.2, 0.2), Point_2(0.3, 0.3), Point_2(0.4, 0.4) };

  // Create a lambda function with affine coordinates.
  // This implementation is based on the following paper:
  // S. Waldron. Affine generalized barycentric coordinates.
  // Jaen Journal on Approximation, 3(2):209-226, 2011.
  const auto affine = [](
    const Points_2& vertices, 
    const Point_2& query, 
    Output_iterator coordinates,
    Kernel traits) { 
    
    const std::size_t n = vertices.size();
    const auto lambda = [](const Point_2& p){ return std::make_pair(p, 1.0); };
    const Point_2 b = CGAL::barycenter(
      boost::make_transform_iterator(vertices.begin(), lambda),
      boost::make_transform_iterator(vertices.end()  , lambda), traits);

    MatrixXd V(2, n);
    for (std::size_t i = 0; i < n; ++i) {
      V(0, i) = vertices[i].x() - b.x();
      V(1, i) = vertices[i].y() - b.y();
    }

    const auto A   = V.adjoint();
    const auto mat = V * A;
    const auto inv = mat.inverse();
    
    Point_2 diff; VectorXd vec(2);
    for (std::size_t i = 0; i < n; ++i) {
      const FT x = query.x() - b.x();
      const FT y = query.y() - b.y();
      diff = Point_2(x, y);
      
      vec(0) = V(0, i);
      vec(1) = V(1, i);
      const auto res = inv * vec;

      *(coordinates++) = 
        diff.x() * res(0) + diff.y() * res(1) + 1.0 / double(n);
    }
  };

  // Evaluate affine coordinates for all query points at once.
  std::vector< std::vector<FT> > bs;
  bs.reserve(queries.size());
  CGAL::Barycentric_coordinates::analytic_weights_2(
    vertices, queries, affine, std::back_inserter(bs), traits, point_map);

  // Output affine coordinates.
  std::cout << std::endl << "affine coordinates: " << std::endl << std::endl;
  for (const auto& b : bs) {
    for (std::size_t i = 0; i < b.size() - 1; ++i)
      std::cout << b[i] << ", ";
    std::cout << b[b.size() - 1] << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
