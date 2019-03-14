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

#ifndef CGAL_BARYCENTRIC_COORDINATES_ENUM_2_H
#define CGAL_BARYCENTRIC_COORDINATES_ENUM_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

namespace CGAL {
/*!
  \ingroup PkgBarycentric_coordinates_2
  The namespace `CGAL::Barycentric_coordinates` contains implementations of all 
  generalized barycentric coordinates: 2D, 3D, related enumerations, etc.
*/
namespace Barycentric_coordinates {

/// \name Query Point Location
/// @{

/*!
  Query_point_location defines different possible locations 
  of a query point.
*/
enum class Query_point_location {
    
    /// Location is not known apriori and defined automatically by the algorithm.
    UNSPECIFIED = 0,

    /// Query point is located at the vertex of the polygon.
    ON_VERTEX = 1,

    /// Query point is located on the edge of the polygon.
    ON_EDGE = 2,

    /// Query point is located inside the polygon, excluding the boundary.
    ON_BOUNDED_SIDE = 3,

    /// Query point is located outside the polygon, excluding the boundary.
    ON_UNBOUNDED_SIDE = 4
};

/// @}

/// \name Algorithm Types
/// @{

/*!
  Algorithm_type provides a way to choose a type of the algorithm 
  to compute pointwise barycentric weights.
*/
enum class Algorithm_type {
    
    /// A default slow version of the algorithm, but very precise.
    MAX_PRECISION = 0,

    /// A fast version of the algorithm, which is less precise but much faster.
    MAX_SPEED = 1
};

/// @}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_ENUM_2_H
