namespace CGAL {
namespace Barycentric_coordinates {

/*!
\ingroup PkgBarycentricCoordinates2RefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
discretized domains restricted to a simple polygon.

\cgalHasModel 
- `CGAL::Barycentric_coordinates::Delaunay_domain_2`
*/

class DiscretizedDomain_2 {

public:

  /*!  
    returns the number of vertices in the domain.
  */
  const std::size_t number_of_vertices() const {

  }

  /*!  
    returns a const reference to the vertex with the index `query_index`.
  */
  const Point_2& vertex(
    const std::size_t query_index) const {

  }

  /*!  
    controls if the vertex with the index `query_index` is on the 
    boundary of the polygon.
  */
  const bool is_on_boundary(
    const std::size_t query_index) const {

  }

  /*!  
    fills `neighbors` with the indices of all vertices, which are connected to the 
    vertex with the index `query_index`.
  */
  void operator()(
    const std::size_t query_index, 
    std::vector<std::size_t>& neighbors) {
        
  }

  /*!  
    fills `element` with the indices of all vertices, which are vertices of the 
    element that contains `query`, returns true if `query` belongs to the domain.
  */
  bool locate(
    const Point_2& query, 
    std::vector<std::size_t>& element) {
        
  }
};
}
}
