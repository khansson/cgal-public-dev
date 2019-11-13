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

#ifndef CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL headers.
#include <vector>
#include <utility>
#include <iterator>

// Boost headers.
#include <boost/optional.hpp>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    This function takes the `source` and `target` vertices of a segment and computes 
    the segment coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `coordinates`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

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
    const GeomTraits traits) {

    if (source == target)
      return boost::none;

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto scalar_product_2   = traits.compute_scalar_product_2_object();
    const auto squared_distance_2 = traits.compute_squared_distance_2_object();

    // Project point on the segment.
    const FT opposite_scalar_product = 
    scalar_product_2(query - target, source - target);

    // Compute coordinates.
    const FT b0 = opposite_scalar_product / squared_distance_2(source, target);
    const FT b1 = FT(1) - b0;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;

    return boost::optional<OutputIterator>(coordinates);
  }

  /*!
    This function takes the vertices `p0`, `p1`, and `p2` of a triangle and computes the
    triangle coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `coordinates`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed triangle coordinates.

    \param traits
    An instance of `GeomTraits`.

    \return an optional output iterator.
  */
  template<
  typename OutputIterator,
  typename GeomTraits>
  boost::optional<OutputIterator> triangle_coordinates_2(
    const typename GeomTraits::Point_2& p0, 
    const typename GeomTraits::Point_2& p1, 
    const typename GeomTraits::Point_2& p2, 
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto area_2 = traits.compute_area_2_object();

    if (area_2(p0, p1, p2) == FT(0))
      return boost::none;

    // Compute some related sub-areas.
    const FT A1 = area_2(p1, p2, query);
    const FT A2 = area_2(p2, p0, query);

    // Compute the inverted total area of the triangle.
    const FT inverted_total_area = FT(1) / area_2(p0, p1, p2);

    // Compute coordinates.
    const FT b0 = A1 * inverted_total_area;
    const FT b1 = A2 * inverted_total_area;
    const FT b2 = FT(1) - b0 - b1;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;
    *(coordinates++) = b2;

    return boost::optional<OutputIterator>(coordinates);
  }

  /*!
    This function takes a `query` point and computes boundary barycentric
    coordinates at this point with respect to the vertices of a given `polygon`.
    These coordinates are then returned in `coordinates`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    In all other cases, all the coordinates are set to zero.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `Identity_property_map`.

    \param polygon
    An instance of `Polygon` with vertices of a 2D polygon.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon` 
    to `GeomTraits::Point_2`.

    \return an optional output iterator.
  */
  template<
  typename Polygon,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> boundary_coordinates_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map = VertexMap()) {

    const auto result = 
    internal::locate_wrt_polygon(polygon, query, vertex_map, traits);
    
    const auto location = (*result).first;
    const auto index    = (*result).second;

    return boundary_coordinates_2(
      polygon, query, location, index, coordinates, traits, vertex_map);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename Polygon,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> boundary_coordinates_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const Query_point_location location,
    const std::size_t index,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map = VertexMap()) {

    using FT = typename GeomTraits::FT;

    // Compute coordinates with respect to the query point location.
    switch (location) {
      
      case Query_point_location::ON_VERTEX: {

        CGAL_precondition(index >= 0 && index < polygon.size());
        for (std::size_t i = 0; i < polygon.size(); ++i)
          if (i == index) 
            *(coordinates++) = FT(1);
          else
            *(coordinates++) = FT(0);
        break;
      }
      
      case Query_point_location::ON_EDGE: {  

        CGAL_precondition(index >= 0 && index < polygon.size());
        const std::size_t indexp = (index + 1) % polygon.size();

        const auto& source = get(vertex_map, *(polygon.begin() + index));
        const auto& target = get(vertex_map, *(polygon.begin() + indexp));

        for (std::size_t i = 0; i < polygon.size(); ++i)
          if (i == index)
            segment_coordinates_2(source, target, query, coordinates, traits);
          else
            *(coordinates++) = FT(0);
        break;
      }

      default: {
        internal::get_default(polygon.size(), coordinates);
        break;
      }
    }
    
    return boost::optional<OutputIterator>(coordinates);
  }
  /// \endcond

  /*!
    This function takes a range of `query` points and computes boundary barycentric
    coordinates at each point with respect to the vertices of a given `polygon`. 
    These coordinates are then returned in `coordinates`.

    If a query point from the range does not belong to the polygon's boundary,
    its coordinates are set to zero.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam QueryRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `Identity_property_map`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `Identity_property_map`.

    \param polygon
    An instance of `Polygon` with vertices of a 2D polygon.

    \param queries
    An instance of `QueryRange` with 2D query points.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon` 
    to `GeomTraits::Point_2`.

    \param point_map
    An instance of `PointMap` that maps an item from `queries` 
    to `GeomTraits::Point_2`.

    \return an optional output iterator.
  */
  template<
  typename Polygon,
  typename QueryRange,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2>,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> boundary_coordinates_2(
    const Polygon& polygon,
    const QueryRange& queries,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map = VertexMap(),
    const PointMap point_map = PointMap()) {

    // Typedefs.
    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    // Compute coordinates for all query points.
    std::vector<FT> b;
    b.reserve(polygon.size());

    for (const auto& item : queries) {
      const Point_2& query = get(point_map, item); b.clear();
      boundary_coordinates_2(
        polygon, query,
        std::back_inserter(b),
        traits, vertex_map);
      
      *(coordinates++) = b;
    }
    
    return boost::optional<OutputIterator>(coordinates);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
