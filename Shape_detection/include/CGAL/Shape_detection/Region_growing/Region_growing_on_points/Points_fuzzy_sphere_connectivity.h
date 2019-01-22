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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_FUZZY_SPHERE_CONNECTIVITY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_FUZZY_SPHERE_CONNECTIVITY_H

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
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_maps.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints
    \brief Fuzzy sphere search for neighbors on a set of `Point_2` or `Point_3`.
    \tparam Traits Model of `Kernel`
    \tparam InputRange An arbitrary range with user-defined items, given an IndexToItemMap is provided. 
    The default one is random access.
    \tparam PointMap An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.
    \cgalModels `RegionGrowingConnectivity`
  */
  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap>
  class Points_fuzzy_sphere_connectivity {

  public:
            
    /// \name Types
    /// @{

    using Traits = GeomTraits;

    using Input_range = InputRange;
    ///< An arbitrary range with user-defined items. The default implementation is random access.

    using Point_map = PointMap;
    ///< An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.

    using Point = typename Point_map::value_type;
    ///< Point type, can only be `Point_2` or `Point_3`.

    using FT = typename Traits::FT;
    ///< Number type.

    ///< \cond SKIP_IN_MANUAL
    using Index_to_point_map = 
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value, 
      CGAL::Search_traits_2<Traits>, 
      CGAL::Search_traits_3<Traits> >::type;
                    
    using Search_traits = 
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;
      
    using Splitter = 
    CGAL::Sliding_midpoint<Search_traits>;
      
    using Fuzzy_sphere 
    = CGAL::Fuzzy_sphere<Search_traits>;
      
    using Tree 
    = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;
    ///< \endcond
                
    /// @}

    /// \name Initialization
    /// @{

    /*!
      The constructor that takes a set of items, provided a point_map to access a `Point` from an item, and a search sphere radius.
    */
    Points_fuzzy_sphere_connectivity(
      const Input_range& input_range, 
      const FT search_radius = FT(1), 
      const Point_map point_map = Point_map()) :
    m_input_range(input_range),
    m_search_radius(search_radius),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      m_tree.build();
      CGAL_precondition(search_radius > FT(0));
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      From a query item with the index `query_index`, this function creates a search sphere centered at this item.
      It then uses `CGAL::Kd_tree::search()` to look for the neighbors of the given query and push their indices to `neighbors`.
      \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
    */
    template<typename OutputIterator>
    void get_neighbors(
      const std::size_t query_index, 
      OutputIterator neighbors) const {
                
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());
      
      const Fuzzy_sphere sphere(
        query_index, 
        m_search_radius, 
        FT(0), 
        m_tree.traits());

      m_tree.search(neighbors, sphere);
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
    
    const FT m_search_radius;

    const Point_map m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Tree m_tree;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_FUZZY_SPHERE_CONNECTIVITY_H
