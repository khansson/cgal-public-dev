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

#ifndef CGAL_BARYCENTRIC_COORDINATES_MEAN_VALUE_WEIGHTS_2_H
#define CGAL_BARYCENTRIC_COORDINATES_MEAN_VALUE_WEIGHTS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  template<
  typename VertexRange,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Mean_value_weights_2 {

    public:
      using Vertex_range = VertexRange;
      using Traits = GeomTraits;
      using Vertex_map = VertexMap;
      
      using Point_2 = typename GeomTraits::Point_2;

      Mean_value_weights_2() {
        
      }

      template<typename OutputIterator>
      boost::optional<OutputIterator> operator()(
        const Point_2& p, 
        OutputIterator weights) const {

      }

      bool is_valid_point(const Point_2& query) const {

      }

      bool is_boundary_point(const Point_2& query) const {

      }

    private:
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_MEAN_VALUE_WEIGHTS_2_H
