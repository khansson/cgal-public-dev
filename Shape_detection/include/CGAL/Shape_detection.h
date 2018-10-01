// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// Author(s) : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

/**
* \ingroup PkgShapeDetection
* \file CGAL/Shape_detection.h
* A header file that includes all headers of this package.
*/

#ifndef CGAL_SHAPE_DETECTION_H
#define CGAL_SHAPE_DETECTION_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/RANSAC_on_points_3/Efficient_RANSAC_traits.h>
#include <CGAL/Shape_detection/RANSAC_on_points_3/Efficient_RANSAC.h>
#include <CGAL/Shape_detection/RANSAC_on_points_3/Octree.h>
#include <CGAL/Shape_detection/RANSAC_on_points_3/Primitives.h>
#include <CGAL/Shape_detection/RANSAC_on_points_3/Property_maps.h>

#include <CGAL/Region_growing/Region_growing.h>
#include <CGAL/Region_growing/Region_growing_traits.h>
#include <CGAL/Region_growing/Region_growing_on_points.h>
#include <CGAL/Region_growing/Region_growing_on_mesh.h>

// #include <CGAL/Region_growing/Region_growing_on_points/Fuzzy_sphere_connectivity.h>
// #include <CGAL/Region_growing/Region_growing_on_points/Nearest_neighbor_connectivity.h>
// #include <CGAL/Region_growing/Region_growing_on_points/Propagation_conditions_2.h>
// #include <CGAL/Region_growing/Region_growing_on_points/Propagation_conditions_3.h>

// #include <CGAL/Region_growing/Region_growing_on_mesh/Mesh_connectivity.h>
// #include <CGAL/Region_growing/Region_growing_on_mesh/Propagation_conditions.h>

#include <CGAL/Shape_detection/regularize_planes.h>

#endif // CGAL_SHAPE_DETECTION_H