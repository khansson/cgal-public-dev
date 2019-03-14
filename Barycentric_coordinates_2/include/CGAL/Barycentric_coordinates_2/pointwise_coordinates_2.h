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
    is a model of `Kernel`.

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
    This function takes the three vertices of a triangle and computes the
    triangle coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `coordinates`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `Kernel`.

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
  typename OutputIterator,
  typename GeomTraits>
  boost::optional<OutputIterator> triangle_coordinates_2(
    const typename GeomTraits::Point_2& p1, 
    const typename GeomTraits::Point_2& p2, 
    const typename GeomTraits::Point_2& p3, 
    const typename GeomTraits::Point_2& query,
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
    This function takes a query point and computes the boundary barycentric
    coordinates at this point with respect to the vertices of a given polygon. 
    These coordinates are then returned in `coordinates`.

    These coordinates can be computed if and only if query point is at the vertex 
    of the polygon or belongs to one of its edges, otherwise all coordinates are
    set to zero.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \param polygon
    An instance of `Polygon` with vertices of a 2D polygon.

    \param query
    A query point.

    \param location
    A query point location with respect to the polygon.

    \param index
    An index of the polygon vertex or edge to which the 
    query point belongs.

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
    const Query_point_location location,
    const std::size_t index,
    OutputIterator coordinates,
    const GeomTraits traits = GeomTraits(),
    const VertexMap vertex_map = VertexMap()) {

    // Compute coordinates with respect to the query point location.
    switch (location) {
      
      case Query_point_location::VERTEX: {
        for (std::size_t i = 0; i < polygon.size(); ++i)
          if (i == index) 
            *(coordinates++) = FT(1);
          else
            *(coordinates++) = FT(0);
        break;
      }
      
      case Query_point_location::EDGE: {  
        const std::size_t indexp = (index + 1) % polygon.size();

        const Point_2& source = get(vertex_map, *(polygon.begin() + index));
        const Point_2& target = get(vertex_map, *(polygon.begin() + indexp));

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

  /*!
    This function takes a range of query points and computes the boundary barycentric
    coordinates at each point with respect to the vertices of a given polygon. 
    These coordinates are then returned in `coordinates`.

    This function is using internally the function above.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam QueryRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

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
    const GeomTraits traits = GeomTraits(),
    const VertexMap vertex_map = VertexMap(),
    const PointMap point_map = PointMap()) {

    // Typedefs.
    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    // Compute coordinates for all query points.
    std::vector<FT> b;
    b.reserve(polygon.size());

    for (auto it = queries.begin(); it != queries.end(); ++it) {
      const Point_2& query = get(point_map, *it); b.clear();

      // Get boundary point location and compute coordinates.
      const auto result = 
      internal::is_boundary_point(polygon, query, vertex_map);
      
      boundary_coordinates_2(
        polygon, query,
        (*result).first, 
        (*result).second,
        std::back_inserter(b),
        traits, vertex_map);
      
      *(coordinates++) = b;
    }
    return boost::optional<OutputIterator>(coordinates);
  }

  /*!
    This function takes a range of query points and computes the chosen 
    barycentric weights at each point with respect to the vertices of a given polygon. 
    These weights are then returned in `output`.

    These weights can be computed if and only if query point is not on the 
    polygon's boundary, otherwise all weights are set to zero.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam QueryRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam Weights
    is a model of `Pointwise_weights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \param polygon
    An instance of `Polygon` with vertices of a 2D polygon.

    \param queries
    An instance of `QueryRange` with 2D query points.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

    \param output 
    An output iterator that stores the computed weights.

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
  typename Weights,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2>,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> pointwise_weights_2(
    const Polygon& polygon,
    const QueryRange& queries,
    const Weights& weights,
    OutputIterator output,
    const GeomTraits traits = GeomTraits(),
    const VertexMap vertex_map = VertexMap(),
    const PointMap point_map = PointMap()) {

    // Typedefs.
    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    // Compute weights for all query points.
    std::vector<FT> w;
    w.reserve(polygon.size());

    for (auto it = queries.begin(); it != queries.end(); ++it) {
      const Point_2& query = get(point_map, *it); w.clear();

      // If point is not valid, skip it.
      if (!weights.is_valid_point(query)) {
        internal::get_default(polygon.size(), std::back_inserter(w));
        *(output++) = w;
        continue;
      }

      // If point is on the boundary, skip it.
      const auto result = weights.is_boundary_point(query);
      if (result) 
        internal::get_default(polygon.size(), std::back_inserter(w));
      else {
        weights(
          query, 
          std::back_inserter(w));
      }
      *(output++) = w;
    }
    return boost::optional<OutputIterator>(output);
  }

  /*!
    This function takes a range of query points and computes the chosen 
    barycentric weights at each point with respect to the vertices of a given polygon. 
    These weights are then normalized and returned in `coordinates`.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam QueryRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam Weights
    is a model of `Pointwise_weights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

    \param polygon
    An instance of `Polygon` with vertices of a 2D polygon.

    \param queries
    An instance of `QueryRange` with 2D query points.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

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
  typename Weights,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2>,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> pointwise_coordinates_2(
    const Polygon& polygon,
    const QueryRange& queries,
    const Weights& weights,
    OutputIterator coordinates,
    const GeomTraits traits = GeomTraits(),
    const VertexMap vertex_map = VertexMap(),
    const PointMap point_map = PointMap()) {

    // Typedefs.
    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    // Compute coordinates for all query points.
    std::vector<FT> b;
    b.reserve(polygon.size());

    for (auto it = queries.begin(); it != queries.end(); ++it) {
      const Point_2& query = get(point_map, *it); b.clear();

      // If point is not valid, skip it.
      if (!weights.is_valid_point(query)) {
        internal::get_default(polygon.size(), std::back_inserter(b));
        *(coordinates++) = b;
        continue;
      }

      // Otherwise compute either boundary or normalized coordinates.
      const auto result = weights.is_boundary_point(query);
      if (result) { 
        boundary_coordinates_2(
          polygon, query,
          (*result).first, 
          (*result).second,
          std::back_inserter(b),
          traits, vertex_map);
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

  /*!
    This function takes a query point and computes the chosen barycentric
    weights at this point with respect to the vertices of a given polygon. 
    These weights are then normalized and returned in `coordinates`.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam Weights
    is a model of `Pointwise_weights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`. %Default is `CGAL::Identity_property_map`.

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
  typename Weights,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional<OutputIterator> pointwise_coordinates_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const Weights& weights,
    OutputIterator coordinates,
    const GeomTraits traits = GeomTraits(),
    const VertexMap vertex_map = VertexMap()) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // If point is not valid, skip it.
    if (!weights.is_valid_point(query)) {
      internal::get_default(polygon.size(), coordinates);
      return boost::optional<OutputIterator>(coordinates);
    }
    
    // Otherwise compute either boundary or normalized coordinates.
    const auto result = weights.is_boundary_point(query);
    if (result) {
      boundary_coordinates_2(
        polygon, query,
        (*result).first, 
        (*result).second,
        coordinates,
        traits, vertex_map);
    } else {
      std::vector<FT> b;
      b.reserve(polygon.size());

      weights(
        query, 
        std::back_inserter(b));
      
      internal::normalize(b);
      for (const auto& value : b)
        *(coordinates++) = value;
    }
    return boost::optional<OutputIterator>(coordinates);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_POINTWISE_COORDINATES_2_H
