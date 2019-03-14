/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
barycentric weight functions that can be computed per point.

\cgalHasModel 
- `CGAL::Barycentric_coordinates::Wachspress_weights_2`,
- `CGAL::Barycentric_coordinates::Mean_value_weights_2`,
- `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2`,
- `CGAL::Barycentric_coordinates::Harmonic_coordinates_2`
*/

class PointwiseWeights2 {

public:

  /*!  
    fills `weights` with the corresponding barycentric weights 
    computed at the `query` point. 

    \return optional output iterator
  */
  template<typename OutputIterator>
  boost::optional<OutputIterator> operator()(
    const Point_2& query, 
    OutputIterator weights) {
        
  }

  /*!  
    checks if the corresponding barycentric weights are well-defined 
    at the `query` point.
  */
  bool is_valid_point(const Point_2& query) {

  }

  /*!  
    checks if the `query` point is on the polygon's boundary.
  */
  bool is_boundary_point(const Point_2& query) {

  }
};
