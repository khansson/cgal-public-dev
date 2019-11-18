// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
  typename OutputIterator,
  typename GeomTraits>
  boost::optional<OutputIterator> segment_coordinates_2(
    const typename GeomTraits::Point_2& source, 
    const typename GeomTraits::Point_2& target, 
    const typename GeomTraits::Point_2& query, 
    OutputIterator coordinates,
    const GeomTraits traits) {

    return internal::linear_coordinates_2(
      source, target, query, coordinates, traits);
  }

  /*!
    This function takes the `source` and `target` vertices of a segment and computes 
    the segment coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `coordinates`.

    This function infers a traits class from the `Point_2` class and calls the
    generic function above.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param source
    The source vertex of a segment.

    \param target
    The target vertex of a segment.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed segment coordinates.

    \return an optional output iterator.
  */
  template<
  typename Point_2, 
  typename OutputIterator>
  boost::optional<OutputIterator> segment_coordinates_2(
    const Point_2& source, 
    const Point_2& target, 
    const Point_2& query, 
    OutputIterator coordinates) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return segment_coordinates_2(
      source, target, query, coordinates, GeomTraits());
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

    return coordinates;
  }

  /*!
    This function takes the vertices `p0`, `p1`, and `p2` of a triangle and computes the
    triangle coordinates at a given `query` point with respect to these vertices.
    The coordinates are returned in `coordinates`.

    This function infers a traits class from the `Point_2` class and calls the
    generic function above.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

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

    \return an optional output iterator.
  */
  template<
  typename Point_2, 
  typename OutputIterator>
  boost::optional<OutputIterator> triangle_coordinates_2(
    const Point_2& p0, 
    const Point_2& p1, 
    const Point_2& p2, 
    const Point_2& query,
    OutputIterator coordinates) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return triangle_coordinates_2(
      p0, p1, p2, query, coordinates, GeomTraits());
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
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \param polygon
    An instance of `Polygon` with vertices of a simple polygon.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon` 
    to `GeomTraits::Point_2`.

    \return an optional output iterator with the tag indicating if `query` 
    belongs to the polygon's boundary.
  */
  template<
  typename Polygon,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap>
  boost::optional< std::pair<OutputIterator, bool> > boundary_coordinates_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using Point_2 = typename GeomTraits::Point_2;
    std::vector<Point_2> poly;
    poly.reserve(polygon.size());
    for (const auto& item : polygon)
      poly.push_back(get(vertex_map, item));

    const auto result = 
    internal::locate_wrt_polygon_2(poly, query, traits);
    const auto location = (*result).first;
    const auto index    = (*result).second;

    return internal::boundary_coordinates_2(
      poly, query, location, index, coordinates, traits);
  }

  /*!
    This function takes a `query` point and computes boundary barycentric
    coordinates at this point with respect to the vertices of a given `polygon`.
    These coordinates are then returned in `coordinates`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    In all other cases, all the coordinates are set to zero.

    This function infers a traits class from the `Point_2` class and calls the
    generic function above.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`.

    \param polygon
    An instance of `Polygon` with vertices of a simple polygon.

    \param query
    A query point.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon` to `Point_2`.

    \return an optional output iterator with the tag indicating if `query` 
    belongs to the polygon's boundary.
  */
  template<
  typename Polygon,
  typename Point_2,
  typename OutputIterator,
  typename VertexMap>
  boost::optional< std::pair<OutputIterator, bool> > boundary_coordinates_2(
    const Polygon& polygon,
    const Point_2& query,
    OutputIterator coordinates,
    const VertexMap vertex_map) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return boundary_coordinates_2(
      polygon, query, coordinates, GeomTraits(), vertex_map);
  }

  /*!
    This function takes a range of `query` points and computes boundary barycentric
    coordinates at each point with respect to the vertices of a given `polygon`. 
    These coordinates are then returned in `coordinates`.

    If a query point from the range does not belong to the polygon's boundary,
    its coordinates are set to zero.

    This function calls the function above on a range of query points.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam QueryRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`.

    \param polygon
    An instance of `Polygon` with vertices of a simple polygon.

    \param queries
    An instance of `QueryRange` with query points.

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
  typename VertexMap,
  typename PointMap>
  boost::optional<OutputIterator> boundary_coordinates_2(
    const Polygon& polygon,
    const QueryRange& queries,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map,
    const PointMap point_map) {

    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    std::vector<Point_2> poly;
    poly.reserve(polygon.size());
    for (const auto& item : polygon)
      poly.push_back(get(vertex_map, item));

    std::vector<FT> b;
    b.reserve(polygon.size());

    for (const auto& item : queries) {
      const auto& query = get(point_map, item);

      const auto result = 
      internal::locate_wrt_polygon_2(poly, query, traits);
      const auto location = (*result).first;
      const auto index    = (*result).second;

      b.clear();
      internal::boundary_coordinates_2(
        poly, query, location, index, std::back_inserter(b), traits);
      *(coordinates++) = b;
    }
    return coordinates;
  }

  /*!
    This function takes a range of `query` points and computes the chosen 
    barycentric `weights` at each point with respect to the given `vertices`. 
    These weights are then returned in `output`.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam QueryRange
    is a model of `ConstRange`.

    \tparam Weights
    is a model of `AnalyticWeights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`.

    \param vertices
    An instance of `VertexRange` with vertices.

    \param queries
    An instance of `QueryRange` with query points.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

    \param output 
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param point_map
    An instance of `PointMap` that maps an item from `queries` 
    to `GeomTraits::Point_2`.

    \return an optional output iterator.
  */
  template<
  typename VertexRange,
  typename QueryRange,
  typename Weights,
  typename OutputIterator,
  typename GeomTraits,
  typename PointMap>
  boost::optional<OutputIterator> analytic_weights_2(
    const VertexRange& vertices,
    const QueryRange& queries,
    Weights& weights,
    OutputIterator output,
    const GeomTraits traits,
    const PointMap point_map) {

    using FT = typename GeomTraits::FT;

    std::vector<FT> w;
    w.reserve(vertices.size());

    for (const auto& item : queries) {
      const auto& query = get(point_map, item);
      w.clear(); weights(
        vertices, query, std::back_inserter(w), traits); 
      *(output++) = w;
    }
    return output;
  }

  /*!
    This function takes a `query` point and computes the chosen barycentric
    `weights` at this point with respect to the given `vertices`. These weights 
    are then normalized and returned in `coordinates`.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam Weights
    is a model of `AnalyticWeights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \param vertices
    An instance of `VertexRange` with vertices.

    \param query
    A query point.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \return an optional output iterator.
  */
  template<
  typename VertexRange,
  typename Point_2,
  typename Weights,
  typename OutputIterator,
  typename GeomTraits>
  boost::optional<OutputIterator> analytic_coordinates_2(
    const VertexRange& vertices,
    const Point_2& query,
    Weights& weights,
    OutputIterator coordinates,
    const GeomTraits traits) {

    using FT = typename GeomTraits::FT;

    std::vector<FT> b;
    b.reserve(vertices.size());

    weights(vertices, query, std::back_inserter(b), traits);
    internal::normalize(b);
    for (const auto& value : b)
      *(coordinates++) = value;
    return coordinates;
  }

  /*!
    This function takes a `query` point and computes the chosen barycentric
    `weights` at this point with respect to the given `vertices`. These weights 
    are then normalized and returned in `coordinates`.

    This function infers a traits class from the `Point_2` class and calls the
    generic function above.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam Weights
    is a model of `AnalyticWeights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param vertices
    An instance of `VertexRange` with vertices.

    \param query
    A query point.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \return an optional output iterator.
  */
  template<
  typename VertexRange,
  typename Point_2,
  typename Weights,
  typename OutputIterator>
  boost::optional<OutputIterator> analytic_coordinates_2(
    const VertexRange& vertices,
    const Point_2& query,
    Weights& weights,
    OutputIterator coordinates) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return analytic_coordinates_2(
      vertices, query, weights, coordinates, GeomTraits());
  }

  /*!
    This function takes a range of `query` points and computes the chosen 
    barycentric `weights` at each point with respect to the given `vertices`. 
    These weights are then normalized and returned in `coordinates`.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam QueryRange
    is a model of `ConstRange`.

    \tparam Weights
    is a model of `AnalyticWeights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `std::vector<GeomTraits::FT>`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam PointMap
    is an `LvaluePropertyMap` whose key type is `QueryRange::value_type` and
    value type is `GeomTraits::Point_2`.

    \param vertices
    An instance of `VertexRange` with vertices.

    \param queries
    An instance of `QueryRange` with query points.

    \param weights
    An instance of `Weights` that computes the corresponding 
    barycentric weights.

    \param coordinates 
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param point_map
    An instance of `PointMap` that maps an item from `queries` 
    to `GeomTraits::Point_2`.

    \return an optional output iterator.
  */
  template<
  typename VertexRange,
  typename QueryRange,
  typename Weights,
  typename OutputIterator,
  typename GeomTraits,
  typename PointMap>
  boost::optional<OutputIterator> analytic_coordinates_2(
    const VertexRange& vertices,
    const QueryRange& queries,
    Weights& weights,
    OutputIterator coordinates,
    const GeomTraits traits,
    const PointMap point_map) {

    using FT = typename GeomTraits::FT;

    std::vector<FT> b;
    b.reserve(vertices.size());

    for (const auto& item : queries) {
      const auto& query = get(point_map, item);
      b.clear(); weights(vertices, query, std::back_inserter(b), traits);
      internal::normalize(b); *(coordinates++) = b;
    }
    return coordinates;
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
