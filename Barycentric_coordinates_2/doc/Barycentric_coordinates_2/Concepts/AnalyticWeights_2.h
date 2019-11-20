/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
barycentric weight functions, which can be computed or evaluated analytically.

\cgalHasModel 
- `CGAL::Barycentric_coordinates::Wachspress_weights_2`
- `CGAL::Barycentric_coordinates::Mean_value_weights_2`
- `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2`
- `CGAL::Barycentric_coordinates::Harmonic_coordinates_2`
*/

class AnalyticWeights_2 {

public:

  /*!  
    fills `weights` with the barycentric weights computed at the `query` point 
    with respect to the given `vertices`, whereas all geometric predicates and 
    constructions are defined in `traits`.
  */
  template<
  typename VertexRange,
  typename Point_2,
  typename OutputIterator,
  typename GeomTraits>
  boost::optional<OutputIterator> operator()(
    const VertexRange& vertices,
    const Point_2& query, 
    OutputIterator weights,
    GeomTraits traits) {
        
  }
};
