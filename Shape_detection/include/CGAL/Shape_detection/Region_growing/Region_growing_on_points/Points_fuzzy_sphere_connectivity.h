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

// #include <CGAL/license/Shape_detection.h>

// STL includes.
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

    \brief Kd tree based fuzzy sphere search for neighbors on a set of `Point_2` or `Point_3`.

    This class uses a Kd tree to search for all points, which belong to a sphere
    of the fixed radius centered at the query point, and thus being its 
    direct neighbors.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam InputRange 
    is a model of `ConstRange`. Its iterator type is `RandomAccessIterator`. 
    Its value type depends on the item type used in Region Growing, 
    for example it can be `CGAL::Point_2`, `CGAL::Point_3`, or 
    or any user-defined type.

    \tparam PointMap 
    is an `LvaluePropertyMap` that maps to `CGAL::Point_2` or `CGAL::Point_3`.

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

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point = typename Point_map::value_type;
    /// \endcond

    /// Number type.
    using FT = typename GeomTraits::FT;

    /// \cond SKIP_IN_MANUAL
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
    /// \endcond
                
    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param input_range 
      An instance of an `InputRange` container with 2D or 3D points.

      \param search_radius 
      Fixed radius of the fuzzy sphere used for searching neighbors.

      \param point_map
      An instance of an `LvaluePropertyMap` that maps an item from `input_range` 
      to `CGAL::Point_2` or to `CGAL::Point_3`.
    */
    Points_fuzzy_sphere_connectivity(
      const InputRange& input_range, 
      const FT search_radius = FT(1), 
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_search_radius(search_radius),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      CGAL_precondition(m_input_range.size() > 0);

      m_tree.build();
      CGAL_precondition(m_search_radius >= FT(0));
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief Returns all points, which are neighbors of a query point.

      This function returns indices of all points, which belong to a sphere
      of the fixed radius `search_radius` centered at the query point with
      the index `query_index`, and thus being its direct neighbors. These
      neighbors are returned in `neighbors`.

      \param query_index
      Index of the query point.

      \param neighbors
      An `std::vector<std::size_t>` with the indices of points, which are 
      neighbors of the point with the index `query_index`.

      Implements the function `RegionGrowingConnectivity::neighbors()`.

      \pre `query_index >= 0 && query_index < total_number_of_points`
    */
    void neighbors(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
                
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());
      
      const Fuzzy_sphere sphere(
        query_index, 
        m_search_radius, 
        FT(0), 
        m_tree.traits());

      m_tree.search(std::back_inserter(neighbors), sphere);
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      Clears all internal data structures.
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
