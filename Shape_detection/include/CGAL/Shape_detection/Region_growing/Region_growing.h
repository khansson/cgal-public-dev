// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_H

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/Iterator_range.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRef
    \brief Region growing algorithm
    \tparam InputRange An arbitrary range with user-defined items.
    \tparam Connectivity_ Model of `RegionGrowingConnectivity`
    \tparam Conditions_ Model of `RegionGrowingPropagationConditions`
  */
  template<class InputRange, class Connectivity_, class Conditions_>
  class Region_growing {

  public:

    /// \name Types
    /// @{

    using Input_range = InputRange;
    ///< An arbitrary range with user-defined items. The range must implement the function size().

    using Connectivity = Connectivity_;

    using Conditions = Conditions_;

    ///< \cond SKIP_IN_MANUAL
    using Visited_items = std::vector<bool>;

    using Running_queue = std::queue<std::size_t>;

    using Items = std::vector<std::size_t>;

    using Regions = std::vector<Items>;
    ///< \endcond

    #ifdef DOXYGEN_RUNNING

      using Regions = unspecified_type;
      
      using Items = unspecified_type;

    #endif

    using Region_range = CGAL::Iterator_range<typename Regions::const_iterator>;
    ///< An `Iterator_range` of the iterators in `CGAL::Shape_detection::Region_growing::Regions`.

    using Item_range = CGAL::Iterator_range<typename Items::const_iterator>;
    ///< An `Iterator_range` of the iterators in `CGAL::Shape_detection::Region_growing::Items`.

    /// @}

    /// \name Initialization
    /// @{

    /*!
      The constructor requires an input range and instances 
      of the Connectivity class and Conditions class.
    */
    Region_growing(
      const Input_range& input_range, 
      Connectivity& connectivity, 
      Conditions& conditions) :
    m_input_range(input_range),
    m_connectivity(connectivity),
    m_conditions(conditions)
    { }

    /// @}

    /// \name Detection 
    /// @{

    /*!
      Perform the region growing algorithm over the input range, 
      using the Connectivity class to find neighbors 
      and the Conditions class to validate regions.
    */
    void detect() {

      clear();
      Items region;

      for (std::size_t seed_index = 0; 
        seed_index < m_input_range.size(); 
        ++seed_index) {
                    
        // Try to grow a new region from the index of the seed item.
        if (!m_visited[seed_index]) {
          propagate(seed_index, region);

          // Check global conditions.
          if (!m_conditions.is_valid_region(region)) 
            revert(region);
          else 
            m_regions.push_back(region);
        }
      }
      m_output_regions = 
      Region_range(m_regions.begin(), m_regions.end());

      // Return indices of all unassigned items.
      for (std::size_t item_index = 0; 
        item_index < m_input_range.size(); 
        ++item_index) {
                    
        if (!m_visited[item_index]) 
          m_unassigned.push_back(item_index);
      }

      m_output_unassigned 
      = Item_range(m_unassigned.begin(), m_unassigned.end());
    }

    /// @}

    /// \name Access
    /// @{  

    /*!
      Return a pair of begin iterator and pass-the-end iterator of the container with found regions. 
      If the function `CGAL::Shape_detection::Region_growing::detect()` has not been called, 
      the first and second of the pair will be the same, which implies an empty container.
    */
    const Region_range& regions() const {
      return m_output_regions;
    }

    /*!
      Return a pair of begin iterator and pass-the-end iterator of the container with indices of all unassigned items. 
      If the function `CGAL::Shape_detection::Region_growing::detect()` has not been called, 
      the first and second of the pair will be the same, which implies an empty container.
    */
    const Item_range& unassigned_items() const {
      return m_output_unassigned;
    }

    /*!
      Return the number of found regions.
    */
    const std::size_t number_of_regions() const {
      return m_regions.size();
    }

    /*!
      Return the number of unassigned items.
    */
    const std::size_t number_of_unassigned_items() const {
      return m_unassigned.size();
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      Clear all internal data structures.
    */
    void clear() {
                
      m_visited.clear();
      m_regions.clear();

      m_unassigned.clear();
      m_visited.resize(m_input_range.size(), false);
    }

    /// @}

  private:

    void propagate(const std::size_t seed_index, Items& region) {
      region.clear();

      // Use two queues, while running on this queue, push to the other queue;
      // When the queue is done, update the shape of the current region and swap to the other queue;
      // depth_index is the index of the queue we are using.
      Running_queue running_queue[2];
      bool depth_index = 0;

      // Once the index of an item is pushed to the queue, it is pushed to the region too.
      m_visited[seed_index] = true;
      running_queue[depth_index].push(seed_index);
      region.push_back(seed_index);

      // Update internal properties of the propagating region.
      m_conditions.update(region);

      while (
        !running_queue[depth_index].empty() || 
        !running_queue[!depth_index].empty()) {

        // Call the next item index of the queue and remove it from the queue.
        const std::size_t item_index = running_queue[depth_index].front();
        running_queue[depth_index].pop();

        // Get neighbors of the current item.
        Items neighbors;
        m_connectivity.get_neighbors(item_index, std::back_inserter(neighbors));

        // Visit the neighbors.
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
          const std::size_t neighbor_index = neighbors[i];
                        
          if (
            !m_visited[neighbor_index] && 
            m_conditions.belongs_to_region(neighbor_index, region)) {

            // Add this neighbor to the other queue so that we can visit it later.
            m_visited[neighbor_index] = true;
            running_queue[!depth_index].push(neighbor_index);
            region.push_back(neighbor_index);
          }
        }

        // Update internal properties of the propagating region.
        if (running_queue[depth_index].empty()) {

          m_conditions.update(region);
          depth_index = !depth_index;
        }
      }
    }

    void revert(const Items& region) {
      for (std::size_t i = 0; i < region.size(); ++i)
        m_visited[region[i]] = false;
    }

    // Fields.
    const Input_range& m_input_range;
    Connectivity& m_connectivity;
    Conditions& m_conditions;

    Visited_items m_visited;
    Regions m_regions;
    Items m_unassigned;

    Region_range m_output_regions = 
    Region_range(m_regions.begin(), m_regions.end());

    Item_range m_output_unassigned = 
    Item_range(m_unassigned.begin(), m_unassigned.end());
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
