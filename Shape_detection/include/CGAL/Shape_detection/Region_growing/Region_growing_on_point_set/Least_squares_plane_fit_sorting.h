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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utilities.h>
#include <CGAL/Shape_detection/Region_growing/internal/property_maps.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*! 
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Best least squares plane fit sorting of 3D points.

    This class allows to sort indices of 3D points, where the sorting is
    based on the quality of the local best least squares plane fit. 
    The plane is fitted to each point and its neighbors found via the 
    `Connectivity` class. The points are sorted in the decreasing quality order,
    that is the first index - the best quality.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam InputRange 
    is a model of `ConstRange`. Its iterator type is `RandomAccessIterator`. 
    Its value type depends on the item type used in Region Growing, 
    for example it can be `std::pair<CGAL::Point_3, int>` 
    or any user-defined type.

    \tparam Connectivity 
    is a model of `RegionGrowingConnectivity`.

    \tparam PointMap 
    is an `LvaluePropertyMap` that maps to `CGAL::Point_3`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery,
  typename PointMap>
  class Least_squares_plane_fit_sorting {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Point_map = PointMap;
    using Seed_map = internal::Seed_property_map;
    /// \endcond
    
    #ifdef DOXYGEN_RUNNING
      /// Property map that returns the quality ordered seed indices of the points.
      typedef unspecified_type Seed_map;
    #endif

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param input_range 
      An instance of an `InputRange` container with 3D points.

      \param connectivity 
      An instance of the `Connectivity` class that is used
      internally to access point neighbors.

      \param point_map
      An instance of an `LvaluePropertyMap` that maps an item from `input_range` 
      to `CGAL::Point_3`.

      \pre `input_range.size() > 0`
    */
    Least_squares_plane_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const PointMap point_map = Point_map()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(point_map) { 
      
      CGAL_precondition(input_range.size() > 0);
      
      m_order.resize(m_input_range.size());
      for (std::size_t i = 0; i < m_input_range.size(); ++i) 
        m_order[i] = i;
      m_scores.resize(m_input_range.size());
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief Sorts point indices.
    */
    void sort() {
      
      compute_scores();
      CGAL_postcondition(m_scores.size() > 0);

      Compare_scores cmp(m_scores);
      std::sort(m_order.begin(), m_order.end(), cmp);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief Returns the `Seed_map` that allows to access 
      the ordered point indices.
    */
    Seed_map seed_map() {
      return Seed_map(m_order);
    }

    /// @}

  private:

    // Types.
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;
    using Compare_scores = internal::Compare_scores<Local_FT>;

    // Functions.
    void compute_scores() {

      std::vector<std::size_t> neighbors;
      std::vector<Local_point_3> points;

      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        neighbors.push_back(i);

        points.clear();
        for (std::size_t j = 0; j < neighbors.size(); ++j) {

          CGAL_precondition(neighbors[j] >= 0);
          CGAL_precondition(neighbors[j] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + neighbors[j]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == neighbors.size());

        Local_plane_3 fitted_plane;
        Local_point_3 fitted_centroid;

        #ifndef CGAL_EIGEN3_ENABLED
          m_scores[i] = linear_least_squares_fitting_3(
            points.begin(), points.end(), 
            fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 3>());
        #else 
          m_scores[i] = linear_least_squares_fitting_3(
            points.begin(), points.end(), 
            fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 3>());
        #endif
      }
    }

    // Fields.
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    const Point_map m_point_map;
    
    std::vector<std::size_t> m_order;
    std::vector<Local_FT> m_scores;

    const To_local_converter m_to_local_converter;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H
