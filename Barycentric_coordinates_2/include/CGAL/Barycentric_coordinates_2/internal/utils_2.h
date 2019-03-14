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

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

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
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  boost::optional< std::pair<Query_point_location, std::size_t> >
  is_boundary_point(
    const Polygon& polygon, 
    const typename GeomTraits::Point_2& query,
    const VertexMap vertex_map = VertexMap()) {

    
  }

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_UTILS_2_H
