/*!
\ingroup PkgShapeDetectionConcepts
\cgalConcept

A concept that describes the set of types and methods required by the class `CGAL::Shape_detection::Region_growing`.

\cgalHasModel `CGAL::Shape_detection::Points_fuzzy_sphere_connectivity`, `CGAL::Shape_detection::Points_k_nearest_neighbors_connectivity`, and `CGAL::Shape_detection::Polygon_mesh_adjacent_faces_connectivity`
*/

class RegionGrowingConnectivity {

public:
    
  /// Find all items, which are connected to an item with the index `query_index`, and push their indices into `neighbors`.
  void get_neighbors(
    std::size_t query_index, 
    OutputIterator &neighbors) {
        
  }
};
