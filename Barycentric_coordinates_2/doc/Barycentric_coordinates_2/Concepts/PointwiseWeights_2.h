/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
barycentric weight functions, which can be computed per point.

\cgalHasModel 
- `CGAL::Barycentric_coordinates::Wachspress_weights_2`,
- `CGAL::Barycentric_coordinates::Mean_value_weights_2`,
- `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2`,
- `CGAL::Barycentric_coordinates::Harmonic_coordinates_2`
*/

class PointwiseWeights_2 {

public:

  /*!  
    fills `weights` with the barycentric weights computed at the `query` point. 

    \return optional output iterator.
  */
  template<typename OutputIterator>
  boost::optional<OutputIterator> operator()(
    const Point_2& query, 
    OutputIterator weights) {
        
  }

  /*!  
    checks if the barycentric weights are well-defined at the `query` point.

    \return boolean `true` or `false`.
  */
  bool is_valid_point(const Point_2& query) {

  }

  /*!  
    checks if the `query` point is on the polygon's boundary.

    \return an optional pair. The first item in the pair is location
    of the `query` point with respect to the polygon. The second item is
    the index of the polygon's vertex or edge if the `query` point belongs to the
    polygon's boundary. It is std::size_t(-1), if it does not.
  */
  boost::optional< std::pair<Query_point_location, std::size_t> > 
  is_boundary_point(const Point_2& query) {

  }
};
