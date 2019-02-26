/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_detection::Region_growing` 
to control if a given set of items forms a region.

\cgalHasModel 
`CGAL::Shape_detection::Points_2_least_squares_line_fit_conditions`, 
`CGAL::Shape_detection::Points_3_least_squares_plane_fit_conditions`, 
`CGAL::Shape_detection::Polygon_mesh_least_squares_plane_fit_conditions`
*/

class RegionType {

public:

  /*!  
    checks if an item with the index `query_index` belongs to the `region` 
    defined by the items with indices in `region`.

    This function is called each time when trying to add a new item to a region.
  */
  bool is_part_of_region(
    const std::size_t query_index, 
    const std::vector<std::size_t>& region) {
        
  }

  /*!  
    checks if the items with indices in region defines a valid region
    
    This function is called each time when propagation is no longer possible 
    for the given seed item. If it is `true`, the region is accepted, otherwise rejected 
    and all its items become again available for the Region Growing.
  */
  bool is_valid_region(
    const std::vector<std::size_t>& region) {
        
  }

  /*!
    This function is called each time when a new seed item is chosen and each time when 
    we expand current region with several new items. In the first case, the number 
    of items in the region is one.
  */
  void update(
    const std::vector<std::size_t>& region) {
    
  }
};
