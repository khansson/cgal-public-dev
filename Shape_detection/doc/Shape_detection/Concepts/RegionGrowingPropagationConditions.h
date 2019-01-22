/*!
\ingroup PkgShapeDetectionConcepts
\cgalConcept

A concept that describes the set of types and methods required by the class `CGAL::Shape_detection::Region_growing`.

\cgalHasModel `CGAL::Shape_detection::Points_2_least_squares_line_fit_conditions`, `CGAL::Shape_detection::Points_3_least_squares_plane_fit_conditions`, and `CGAL::Shape_detection::Polygon_mesh_least_squares_plane_fit_conditions`
*/

class RegionGrowingPropagationConditions {

public:
    
  /// A local condition that checks if an item with the index `query_index` belongs to a `region`.
  bool belongs_to_region(
    const std::size_t query_index, 
    const ItemRange& region) {
        
  }

  /// A global condition that checks the validity of the whole `region`.
  bool are_valid(
    const ItemRange& region) {
        
  }

  /// Update all necessary conditions for the given `region`.
  void update(
    const ItemRange& region) {
    
  }
};
