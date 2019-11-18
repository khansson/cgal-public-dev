#include <vector>
#include <CGAL/barycenter.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.

// General.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Points_2 = std::vector<Point_2>;

// Triangulation.
using FB  = CGAL::Delaunay_mesh_face_base_2<Kernel>;
using VB  = CGAL::Triangulation_vertex_base_2<Kernel>;
using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS>;

using Vertex_handle = typename CDT::Vertex_handle;

// Mesher.
using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
using Mesher   = CGAL::Delaunay_mesher_2<CDT, Criteria>;

// Solver.
using MatrixFT = Eigen::SparseMatrix<FT>;
using Solver   = Eigen::SimplicialLDLT<MatrixFT>;

// Coordinates.
using Harmonic = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<
  Points_2, CDT, Solver, Kernel>;

int main() {  

  // Construct a simple polygon.
  const Points_2 polygon = {
    Point_2(0.03, 0.05), Point_2(0.07, 0.04), Point_2(0.10, 0.04),
    Point_2(0.14, 0.04), Point_2(0.17, 0.07), Point_2(0.19, 0.09),
    Point_2(0.22, 0.11), Point_2(0.25, 0.11), Point_2(0.27, 0.10),
    Point_2(0.30, 0.07), Point_2(0.31, 0.04), Point_2(0.34, 0.03),
    Point_2(0.37, 0.02), Point_2(0.40, 0.03), Point_2(0.42, 0.04),
    Point_2(0.44, 0.07), Point_2(0.45, 0.10), Point_2(0.46, 0.13),
    Point_2(0.46, 0.19), Point_2(0.47, 0.26), Point_2(0.47, 0.31),
    Point_2(0.47, 0.35), Point_2(0.45, 0.37), Point_2(0.41, 0.38),
    Point_2(0.38, 0.37), Point_2(0.35, 0.36), Point_2(0.32, 0.35),
    Point_2(0.30, 0.37), Point_2(0.28, 0.39), Point_2(0.25, 0.40),
    Point_2(0.23, 0.39), Point_2(0.21, 0.37), Point_2(0.21, 0.34),
    Point_2(0.23, 0.32), Point_2(0.24, 0.29), Point_2(0.27, 0.24),
    Point_2(0.29, 0.21), Point_2(0.29, 0.18), Point_2(0.26, 0.16),
    Point_2(0.24, 0.17), Point_2(0.23, 0.19), Point_2(0.24, 0.22),
    Point_2(0.24, 0.25), Point_2(0.21, 0.26), Point_2(0.17, 0.26),
    Point_2(0.12, 0.24), Point_2(0.07, 0.20), Point_2(0.03, 0.15),
    Point_2(0.01, 0.10), Point_2(0.02, 0.07)
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

  // Create seeds for interior faces.
  std::list<Point_2> list_of_seeds;
  list_of_seeds.push_back(Point_2(0.1, 0.1));

  // Refine this triangulation.
  Mesher mesher(cdt);
  mesher.set_seeds(list_of_seeds.begin(), list_of_seeds.end(), true);
  mesher.set_criteria(Criteria(0.01, 0.01));
  mesher.refine_mesh();

  // Instantiate a sparse solver.
  Solver solver;

  // Compute harmonic coordinates at the vertices of the input triangulation.
  Harmonic harmonic(polygon, cdt, solver);
  harmonic.compute();
  
  std::cout << std::endl << 
    "harmonic coordinates computed and evaluated at: " 
  << std::endl;

  // Evaluate harmonic coordinates at the barycenters of the 
  // finite triangulation faces.
  std::vector<FT> coordinates;
  coordinates.reserve(polygon.size());
  for (auto fh = cdt.finite_faces_begin();
  fh != cdt.finite_faces_end(); ++fh) {

    const Point_2 b = CGAL::barycenter(
      fh->vertex(0)->point(), FT(1),
      fh->vertex(1)->point(), FT(1),
      fh->vertex(2)->point(), FT(1));
    
    coordinates.clear();
    CGAL::Barycentric_coordinates::analytic_coordinates_2(
      polygon, b, harmonic, std::back_inserter(coordinates));
    
    std::cout << b << ": ";
    for (std::size_t i = 0; i < coordinates.size() - 1; ++i)
      std::cout << coordinates[i] << ", ";
    std::cout << coordinates[coordinates.size() - 1] << std::endl;
  }
}
