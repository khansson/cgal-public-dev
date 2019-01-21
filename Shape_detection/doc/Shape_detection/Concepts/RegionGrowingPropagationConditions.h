/*!
\ingroup PkgShapeDetectionConcepts
\cgalConcept

A concept that describes the set of types and methods required by the class `CGAL::Shape_detection::Region_growing`.

\cgalHasModel `CGAL::Shape_detection::Propagation_conditions_on_points_2`, `CGAL::Shape_detection::Propagation_conditions_on_points_3`, and `CGAL::Shape_detection::Propagation_conditions_on_face_graph`
*/

class RegionGrowingPropagationConditions {

public:
    
  /// A local condition that checks if an item with the index `query_index` belongs to a `region`, where `Item_index` is any signed integer and `Region` is a random access container with indices of other items.
  bool belongs_to_region(
    const std::size_t query_index, const Region &region) const {
        
  }

  /// A global condition that checks the validity of the whole `region`, where `Region` is a random access container with item indices.
  bool are_valid(
    const Region &region) const {
        
  }

  /// Update all necessary conditions for the given `region` of the type `Region`, which is a random access container with item indices.
  void update(
    const Region &region) {
    
  }
};
