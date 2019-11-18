/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
barycentric weight functions, which can be computed or evaluated analytically.

\cgalHasModel 
- `Barycentric_coordinates::Wachspress_weights_2`,
- `Barycentric_coordinates::Mean_value_weights_2`,
- `Barycentric_coordinates::Discrete_harmonic_weights_2`,
- `Barycentric_coordinates::Harmonic_coordinates_2`
*/

class AnalyticWeights_2 {

public:

  /*!  
    fills `weights` with the barycentric weights computed at the `query` point 
    with respect to the given `vertices`. All geometric predicates and constructions
    are defined in `traits`.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \param vertices
    An instance of `VertexRange` with vertices.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.
    
    \param traits
    An instance of `GeomTraits`.

    \return an optional output iterator.
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
