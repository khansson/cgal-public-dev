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

// #include <CGAL/license/Shape_detection.h>

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/Iterator_range.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRG
    
    \brief Generic Region Growing algorithm.

    This version of Region Growing algorithm allows to grow regions on a set
    of user-defined items given a way to access neighbors of each item via a 
    `Connectivity` parameter class and control if items should form a region
    via a `Conditions` class.
    
    \tparam InputRange 
    is a model of `ConstRange`. Its iterator type is `RandomAccessIterator`. 
    Its value type depends on the item type used in Region Growing, 
    for example it can be `CGAL::Point_2`, `std::pair<CGAL::Point_3, CGAL::Vector_3>`, 
    or any user-defined type.

    \tparam Connectivity 
    is a model of `RegionGrowingConnectivity`.

    \tparam Conditions
    is a model of `RegionGrowingPropagationConditions`
  */
  template<class InputRange, class Connectivity, class Conditions>
  class Region_growing {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Input_range = InputRange;
    using Input_connectivity = Connectivity;
    using Input_conditions = Conditions;

    using Visited_items = std::vector<bool>;
    using Running_queue = std::queue<std::size_t>;    
    /// \endcond

    /// An `std::vector` with indices of all input items provided via `InputRange`.
    using Items = std::vector<std::size_t>;
    
    /// An `std::vector` that stores all found regions, where each region is of type `Region_growing::Items`.
    using Regions = std::vector<Items>;
    
    /// An `Iterator_range` of the iterators in `Region_growing::Regions`.
    using Region_range = CGAL::Iterator_range<typename Regions::const_iterator>;
    
    /// An `Iterator_range` of the iterators in `Region_growing::Items`.
    using Item_range = CGAL::Iterator_range<typename Items::const_iterator>;
    
    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes the Region Growing algorithm.
      
      \param input_range Contains items, which are planned to be used to grow
      regions upon.

      \param connectivity An instance of the `Connectivity` class that is used
      internally to access item neighbors.

      \param conditions An instance of the `Conditions` class that is used
      internally to control if items should form a region.
    */
    Region_growing(
      const InputRange& input_range, 
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
      \brief Runs the Region Growing algorithm.

      Runs the Region Growing algorithm on the input range with items, 
      using a `Connectivity` class to find neighbors and a `Conditions` class 
      to validate regions.
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
      \brief Returns found regions.      

      Returns an `CGAL::Iterator_range` with a bidirectional iterator whose value type
      is `Region_growing::Items`. It is empty if the function `Region_growing::detect()`
      has not been called.
    */
    const Region_range& regions() const {
      return m_output_regions;
    }

    /*!
      \brief Returns indices of all unassigned items.

      Returns an `CGAL::Iterator_range` with a bidirectional iterator whose value type
      is `std::size_t`. It is empty if the function `Region_growing::detect()`
      has not been called.
    */
    const Item_range& unassigned_items() const {
      return m_output_unassigned;
    }

    /*!
      \brief Returns the number of found regions.
    */
    const std::size_t number_of_regions() const {
      return m_regions.size();
    }

    /*!
      \brief Returns the number of unassigned items.
    */
    const std::size_t number_of_unassigned_items() const {
      return m_unassigned.size();
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      Clears all internal data structures.
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
    Input_connectivity& m_connectivity;
    Input_conditions& m_conditions;

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
