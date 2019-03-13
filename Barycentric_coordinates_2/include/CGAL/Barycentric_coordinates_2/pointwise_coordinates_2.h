// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_COORDINATES_POINTWISE_COORDINATES_2_H
#define CGAL_BARYCENTRIC_COORDINATES_POINTWISE_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL headers.
#include <vector>
#include <iterator>

// Boost headers.
#include <boost/optional.hpp>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    This function takes the `source` and `target` vertices of a segment and computes 
    segment coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `output`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \param source
    The source vertex of a segment.

    \param target
    The target vertex of a segment.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed segment coordinates.

    \param traits
    An instance of `GeomTraits`.

    \pre `source != target`

    \return an optional output iterator.
  */
  template<
  typename GeomTraits, 
  typename OutputIterator>
  boost::optional<OutputIterator> segment_coordinates_2(
    const typename GeomTraits::Point_2& source, 
    const typename GeomTraits::Point_2& target, 
    const typename GeomTraits::Point_2& query, 
    OutputIterator coordinates,
    const GeomTraits traits = GeomTraits()) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto scalar_product_2 = traits.compute_scalar_product_2_object();
    const auto squared_distance_2 = traits.compute_squared_distance_2_object();

    // Preconditions.
    CGAL_precondition(source != target);

    // Project point on the segment.
    const FT opposite_scalar_product = 
    scalar_product_2(query - target, source - target);

    // Compute coordinates.
    const FT b1 = opposite_scalar_product / squared_distance_2(source, target);
    const FT b2 = FT(1) - b1;

    // Return coordinates.
    *(coordinates++) = b1;
    *(coordinates++) = b2;

    return boost::optional<OutputIterator>(coordinates);
  }

  /*!
    This function takes the three vertices of a triangle and computes 
    triangle coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `output`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \param p1
    The first vertex of a triangle.

    \param p2
    The second vertex of a triangle.

    \param p3
    The third vertex of a triangle.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed triangle coordinates.

    \param traits
    An instance of `GeomTraits`.

    \pre `area_2(p1, p2, p3) != 0`

    \return an optional output iterator.
  */
  template<
  typename GeomTraits,
  typename OutputIterator>
  boost::optional<OutputIterator> triangle_coordinates_2(
    const typename Traits::Point_2& p1, 
    const typename Traits::Point_2& p2, 
    const typename Traits::Point_2& p3, 
    const typename Traits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits = GeomTraits()) {
    
    // Number type.
    using FT = typename GeomTraits::FT;
    
    // Functions.
    const auto area_2 = traits.compute_area_2_object();

    // Preconditions.
    CGAL_precondition(area_2(p1, p2, p3) != FT(0));

    // Compute some related sub-areas.
    const FT A2 = area_2(p2, p3, query);
    const FT A3 = area_2(p3, p1, query);

    // Compute the inverted total area of the triangle.
    const FT inverted_total_area = FT(1) / area_2(p1, p2, p3);

    // Compute coordinates.
    const FT b1 = A2 * inverted_total_area;
    const FT b2 = A3 * inverted_total_area;
    const FT b3 = FT(1) - b1 - b2;

    // Return coordinates.
    *(coordinates++) = b1;
    *(coordinates++) = b2;
    *(coordinates++) = b3;

    return boost::optional<OutputIterator>(coordinates);
  }

  /*!
    This function takes a range of query points and computes the chosen 
    barycentric weights at each point. These weights are then normalized 
    and returned in `output`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam QueryRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam Weights
    is a model of `Pointwise_weights_2`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \param queries
    An instance of `QueryRange` with 2D query points.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param point_map
    An instance of `PointMap` that maps an item from `queries` 
    to `GeomTraits::Point_2`.

    \return an optional output iterator.
  */
  template<
  typename GeomTraits,
  typename QueryRange,
  typename OutputIterator,
  typename Weights,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> pointwise_coordinates_2(
    const QueryRange& queries,
    const Weights& weights,
    OutputIterator coordinates,
    const PointMap point_map = PointMap()) {

    // Number type.
    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    std::vector<FT> b;
    for (auto it = queries.begin(); it != queries.end(); ++it) {
      const Point_2& query = get(point_map, *it);

      if (!weights.is_valid_point(query))
        continue;

      b.clear();
      if (weights.is_boundary_point(query)) {
        
        boundary_coordinates_2(
          weights.polygon().begin(),
          weights.polygon().end(),
          query,
          std::back_inserter(b));

      } else {
        
        weights(
          query, 
          std::back_inserter(b));
        internal::normalize(b);
        
      }
      *(coordinates++) = b;
    }
    return boost::optional<OutputIterator>(coordinates);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_POINTWISE_COORDINATES_2_H
