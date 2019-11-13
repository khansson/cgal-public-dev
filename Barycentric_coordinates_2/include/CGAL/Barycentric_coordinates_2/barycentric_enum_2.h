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

#ifndef CGAL_BARYCENTRIC_ENUM_2_H
#define CGAL_BARYCENTRIC_ENUM_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

namespace CGAL {

/*!
  \ingroup PkgBarycentric_coordinates_2
  The namespace `Barycentric_coordinates` contains implementations of all 
  generalized barycentric coordinates: 2D, 3D, related enumerations, etc.
*/
namespace Barycentric_coordinates {

/// \name Query Point Location
/// @{

/*!
  `Query_point_location` defines different locations of a query point
  with respect to a polygon.
*/
enum class Query_point_location {

  /// Query point is located at the vertex of the polygon.
  ON_VERTEX = 0,

  /// Query point is located on the edge of the polygon.
  ON_EDGE = 1,

  /// Query point is located in the polygon's interior.
  ON_BOUNDED_SIDE = 2,

  /// Query point is located in the polygon's exterior.
  ON_UNBOUNDED_SIDE = 3,

  /// Location is unspecified. Leads to all coordinates being set to zero.
  UNSPECIFIED = 4
};

/// @}

/// \name Computation Policies
/// @{

/*!
  `Computation_policy` provides a way to choose an asymptotic time complexity 
  of the algorithm.
*/
enum class Computation_policy {

  /*! 
    Computation is very precise but has typically a quadratic time complexity 
    with respect to the number of the polygon's vertices.
  */
  PRECISE_COMPUTATION = 0,

  /*! 
    Computation has typically a linear time complexity with respect to the 
    number of the polygon's vertices, but may be less precise.
  */
  FAST_COMPUTATION = 1
};

/// @}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_2_H
