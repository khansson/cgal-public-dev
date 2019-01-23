/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used to control if items in a 
given set form a region by the `CGAL::Shape_detection::Region_growing` approach.

\cgalHasModel 
`CGAL::Shape_detection::Points_2_least_squares_line_fit_conditions`, 
`CGAL::Shape_detection::Points_3_least_squares_plane_fit_conditions`, 
`CGAL::Shape_detection::Polygon_mesh_least_squares_plane_fit_conditions`
*/

class RegionGrowingPropagationConditions {

public:

  /*!  
    Controls if an item with the index `query_index` belongs to a `region` that
    is defined as `std::vector` with indices of all currently assigned to it items.

    This function is called each time when trying to add a new item to a region.
  */
  bool belongs_to_region(
    const std::size_t query_index, 
    const std::vector<std::size_t>& region) {
        
  }

  /*!  
    Controls if the `region` with currently assigned items given as their indices 
    is globally valid to be accepted as a final region.
    
    This function is called each time when propagation is no longer possible 
    for the given seed item.
  */
  bool is_valid_region(
    const std::vector<std::size_t>& region) {
        
  }

  /*!
    Updates all necessary internal information that goes along with the given `region` 
    that is defined as `std::vector` with indices of all currently assigned to it items.

    This function is called each time when a new seed item is chosen and 
    propagation for this seed has just started and each time right after checking 
    all current query item's neighbors, if they belong to a region, during the propagation.
  */
  void update(
    const std::vector<std::size_t>& region) {
    
  }
};
