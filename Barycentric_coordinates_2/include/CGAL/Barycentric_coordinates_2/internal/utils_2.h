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
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_COORDINATES_UTILS_2_H
#define CGAL_BARYCENTRIC_COORDINATES_UTILS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

  template<typename GeomTraits> 
  class Default_sqrt {
    
  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const { 
      
      CGAL_precondition(value >= FT(0));
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits, 
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) { 
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) { 
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  void normalize(std::vector<FT>& values) {

    FT sum = FT(0);
    for (const auto& value : values)
      sum += value;
    
    CGAL_assertion(sum != FT(0));
    for (auto& value : values)
      value /= sum;
  }

  template<typename OutputIterator>
  void get_default(const std::size_t size, OutputIterator output) {
    for (std::size_t i = 0; i < size; ++i)
      *(output++) = 0;
  }

  template<
  typename Polygon,
  typename VertexMap,
  typename GeomTraits>
  boost::optional< std::pair<Query_point_location, std::size_t> >
  locate_wrt_polygon(
    const Polygon& polygon, 
    const typename GeomTraits::Point_2& query,
    const VertexMap& vertex_map,
    const GeomTraits& traits) {

    
  }

  enum class Polygon_type {
    
    // Concave polygon = non-convex polygon.
    CONCAVE = 0,

    // This is a convex polygon with collinear vertices.
    WEAKLY_CONVEX = 1,

    // This is a convex polygon without collinear vertices.
    STRICTLY_CONVEX = 2
  };

  template<
  typename Polygon,
  typename VertexMap,
  typename GeomTraits>
  Polygon_type polygon_type(
    const Polygon& polygon,
    const VertexMap& vertex_map,
    const GeomTraits& traits) const {
    
    using Collinear_2 = typename GeomTraits::Collinear_2;
    const Collinear_2 collinear_2 = traits.collinear_2_object();

    CGAL_precondition(polygon.size() >= 3);

    std::vector<Point_2> vertices;
    vertices.reserve(polygon.size());
    for (const auto& item : polygon)
      vertices.push_back(get(vertex_map, item));
    CGAL_assertion(vertices.size() == polygon.size());

    // First, test the polygon on convexity.
    if (CGAL::is_convex_2(vertices.begin(), vertices.end(), traits)) {

      // Test all the consequent triplets of the polygon's vertices on collinearity.
      // In case we find at least one, return WEAKLY_CONVEX polygon.
      const std::size_t n = vertices.size();
      for (std::size_t i = 0; i < n; ++i) {
      
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
      
        if (collinear_2(vertices[im], vertices[i], vertices[ip]))
          return Polygon_type::WEAKLY_CONVEX;
      }
      // Otherwise, return STRICTLY_CONVEX polygon.
      return Polygon_type::STRICTLY_CONVEX;
    }
    // Otherwise, return CONCAVE polygon.
    return Polygon_type::CONCAVE;
  }

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_UTILS_2_H
