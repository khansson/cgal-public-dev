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

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <queue>
#include <vector>

// Boost headers.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_maps.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRG
    
    \brief Main class/entry point for running the region growing algorithm.

    This version of the region growing algorithm allows to grow regions on a set
    of user-defined items given a way to access neighbors of each item via a 
    `Connectivity` parameter class and control if items should form a region
    via a `Conditions` class. The `SeedMap` property map allows to define the
    order of items that is which items are used first as seed items to grow
    regions from. It also allows to skip items that should not be used at all.
    
    \tparam InputRange 
    is a model of `ConstRange`. Its iterator type is `RandomAccessIterator`. 

    \tparam Connectivity 
    is a model of `RegionGrowingConnectivity`.

    \tparam Conditions
    is a model of `RegionGrowingPropagationConditions`

    \tparam SeedMap
    is an `LvaluePropertyMap` that maps the `std::size_t` index of an item 
    in `input_range` to an `std::size_t` index of this item in the region growing 
    processing queue. The default one is the identity map.
  */
  template<
  typename InputRange, 
  typename NeighborQuery, 
  typename RegionType,
  typename SeedMap = internal::Identity_seed_property_map>
  class Region_growing {

  public:

    /// \cond SKIP_IN_MANUAL
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Region_type = RegionType;
    using Seed_map = SeedMap;

    using Visited_items = std::vector<bool>;
    using Running_queue = std::queue<std::size_t>;
    using Indices       = std::vector<std::size_t>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes the Region Growing algorithm.
      
      \param input_range 
      Contains items, which are planned to be used to grow regions upon.

      \param neighbor_query 
      An instance of the `Connectivity` class that is used internally to 
      access item neighbors.

      \param region_type 
      An instance of the `Conditions` class that is used internally to 
      control if items should form a region.

      \param seed_map 
      An instance of the `SeedMap` property map that is used internally to 
      set the order of items to grow regions from. If it maps 
      to `std::size_t(-1)`, then this item is skipped.

      \pre `input_range.size() > 0`
    */
    Region_growing(
      const InputRange& input_range, 
      NeighborQuery& neighbor_query, 
      RegionType& region_type,
      const SeedMap seed_map = SeedMap()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_region_type(region_type),
    m_seed_map(seed_map) { 

      CGAL_precondition(input_range.size() > 0);
    }

    /// @}

    /// \name Detection 
    /// @{

    /*!
      \brief Runs the Region Growing algorithm and  fills an output iterator with the results

      Runs the Region Growing algorithm on the input range with items, 
      using a `Connectivity` class to find neighbors and a `Conditions` class 
      to validate regions. The `SeedMap` property map is used to define the
      seeding order of items inside the algorithm.

      This is a useful function if the found regions should be stored
      in the user-defined container outside the class. It helps to avoid
      copying these data from internal to external storage.

      \tparam OutputIterator 
      is an output iterator whose value type is `Item_indices`.

      \param regions
      An output iterator that stores regions represented as `Item_indices`.

      \warning All functions from Access Section will return empty values.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> detect(OutputIterator regions) {

      clear();
      Indices region;

      // Grow regions.
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t seed_index = get(m_seed_map, i);

        // Skip items that user does not want to use.
        if (seed_index == std::size_t(-1))
          continue;

        // Try to grow a new region from the index of the seed item.
        if (!m_visited[seed_index]) {
          propagate(seed_index, region);

          // Check global conditions.
          if (!m_region_type.is_valid_region(region)) 
            revert(region);
          else 
            *(regions++) = region;
        }
      }

      return boost::optional<OutputIterator>(regions);
    }

    /// @}

    /// \name Output
    /// @{  

    /*!
      \brief Returns indices of all unassigned items.

      Returns an `CGAL::Iterator_range` with a bidirectional iterator whose value type
      is `std::size_t`. It is empty if the function `Region_growing::detect()`
      has not been called.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> output_unassigned_items(OutputIterator output) const {
      
      // Return indices of all unassigned items.
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t seed_index = get(m_seed_map, i);

        // Skip items that user does not want to use.
        if (seed_index == std::size_t(-1))
          continue;

        if (!m_visited[seed_index]) 
          *(output++) = seed_index;
      }

      return boost::optional<OutputIterator>(output);
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      Clears all internal data structures.
    */
    void clear() {
                
      m_visited.clear();
      m_visited.resize(m_input_range.size(), false);
    }

    /*!
      Releases all memory that is used internally.
    */
    void release_memory() {

      m_visited.clear();
      m_visited.shrink_to_fit();
    }

    /// @}

  private:

    void propagate(const std::size_t seed_index, Indices& region) {
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

      // Update internal properties of the region.
      m_region_type.update(region);

      Indices neighbors;
      while (
        !running_queue[depth_index].empty() || 
        !running_queue[!depth_index].empty()) {

        // Call the next item index of the queue and remove it from the queue.
        const std::size_t item_index = running_queue[depth_index].front();
        running_queue[depth_index].pop();

        // Get neighbors of the current item.
        neighbors.clear();
        m_neighbor_query(item_index, neighbors);

        // Visit all found neighbors.
        for (const std::size_t neighbor_index : neighbors) {

          // Skip items that user does not want to use.
          if (neighbor_index == std::size_t(-1))
            continue;

          if (!m_visited[neighbor_index] && 
            m_region_type.is_part_of_region(neighbor_index, region)) {

            // Add this neighbor to the other queue so that we can visit it later.
            m_visited[neighbor_index] = true;
            running_queue[!depth_index].push(neighbor_index);
            region.push_back(neighbor_index);
          }
        }

        // Update internal properties of the region.
        if (running_queue[depth_index].empty()) {

          m_region_type.update(region);
          depth_index = !depth_index;
        }
      }
    }

    void revert(const Indices& region) {
      for (const std::size_t item_index : region)
        m_visited[item_index] = false;
    }

    // Fields.
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Region_type& m_region_type;
    const Seed_map m_seed_map;

    Visited_items m_visited;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
