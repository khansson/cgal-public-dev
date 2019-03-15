/*!
\ingroup PkgShapeDetectionRGConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_detection::Region_growing` 
to maintain a region.

A region is represented by a set of indices of the items, which are included in 
this region. These indices are stored in `region`.

Note that you cannot modify `region` since it is constructed inside
`CGAL::Shape_detection::Region_growing`. It is provided only to access its items.

\cgalHasModel 
- `CGAL::Shape_detection::Point_set::Least_squares_line_fit_region`, 
- `CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region`, 
- `CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region`
*/

class RegionType {

public:

  /*!  
    checks if the item with the index `query_index` can be added to the `region`.

    This function is called each time when trying to add a new item to a region.
    If it returns `true`, the `query_index` is pushed to the `region`, otherwise 
    it is ignored.
  */
  bool is_part_of_region(
    const std::size_t query_index, 
    const std::vector<std::size_t>& region) {
        
  }

  /*!  
    checks if the `region` satisfies all necessary conditions.
    
    This function is called at the end of each propagation phase. If it returns `true`, 
    the `region` is accepted, otherwise rejected. If the `region` is rejected,
    all its items are released and available for region growing again.
  */
  bool is_valid_region(
    const std::vector<std::size_t>& region) {
        
  }

  /*!
    allows to update any information that is maintained with the `region`.

    This function is called each time when a new seed item is selected. This item 
    is pushed to the `region` and hence the region's size is one. The function is 
    also called periodically when enlarging the `region`.
  */
  void update(
    const std::vector<std::size_t>& region) {
    
  }
};
