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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H

// STL includes.
#include <vector>
#include <typeinfo>
#include <type_traits>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_maps.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints
    \brief K nearest neighbors (kNN) search on a set of `Point_2` or `Point_3`.
    \tparam Traits Model of `Kernel`
    \tparam InputRange An arbitrary range with user-defined items, given an IndexToItemMap is provided. The default one is random access.
    \tparam PointMap An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.
    \cgalModels `RegionGrowingConnectivity`
  */
  template<class GeomTraits, class InputRange, class PointMap>
  class Points_k_nearest_neighbor_connectivity {

  public:

    /// \name Types
    /// @{

    using Traits = GeomTraits;

    using Input_range = InputRange;
    ///< An arbitrary range with user-defined items.

    using Point_map = PointMap;
    ///< An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.

    using Point = typename Point_map::value_type;
    ///< Point type, can only be `Point_2` or `Point_3`.

    ///< \cond SKIP_IN_MANUAL
    using Index_to_point_map = 
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value, 
      CGAL::Search_traits_2<Traits>, 
      CGAL::Search_traits_3<Traits> >::type;

    using Search_traits = 
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;

    using Distance = 
    CGAL::Distance_adapter<
      std::size_t, 
      Index_to_point_map, 
      CGAL::Euclidean_distance<Search_base> >;

    using Splitter = 
    CGAL::Sliding_midpoint<Search_traits>;

    using Search_tree = 
    CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;

    using Neighbor_search = 
    CGAL::Orthogonal_k_neighbor_search<
      Search_traits, 
      Distance, 
      Splitter, 
      Search_tree>;

    using Tree = 
    typename Neighbor_search::Tree;
    ///< \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      The constructor that takes a set of items, provided a point_map 
      to access a `Point` from an item, 
      and a number of nearest neighbors (the value "k" in "kNN").
    */
    Points_k_nearest_neighbor_connectivity(
      const Input_range& input_range, 
      const std::size_t number_of_neighbors = 12, 
      const Point_map point_map = Point_map()) :
    m_input_range(input_range),
    m_number_of_neighbors(number_of_neighbors),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_distance(m_index_to_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      m_tree.build();
      CGAL_precondition(number_of_neighbors > 0);
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
    This function takes the index `query_index` of a query item and 
    returns indices of the k closest items around it. The result is stored in `neighbors`.
    \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
    */
    template<class OutputIterator>
    void get_neighbors(
      const std::size_t query_index, 
      OutputIterator neighbors) const {
      
      CGAL_precondition(query_index < m_input_range.size());

      Neighbor_search neighbor_search(
        m_tree, 
        get(m_index_to_point_map, query_index), 
        m_number_of_neighbors, 
        0, 
        true, 
        m_distance);
                
      for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
        *(neighbors++) = it->first;
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      Clear all internal data structures.
    */
    void clear() {
      m_tree.clear();
    }

    /// @}

  private:

    // Fields.
    const Input_range& m_input_range;
    
    const std::size_t m_number_of_neighbors;

    const Point_map m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Distance m_distance;
    Tree m_tree;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
