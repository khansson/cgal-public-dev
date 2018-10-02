// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s) : Simon Giraudot
//

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/RANSAC_on_points_3/Primitives/Shape_base.h>

#ifndef CGAL_SHAPE_DETECTION_RANSAC_ON_POINTS_3_PLANE_MAP_H
#define CGAL_SHAPE_DETECTION_RANSAC_ON_POINTS_3_PLANE_MAP_H

namespace CGAL {

namespace Shape_detection {

namespace RANSAC {

	/*!
		\ingroup PkgShapeDetection

		Property map that associates a detected plane object
		(`CGAL::Shape_detection::Plane`) to a `CGAL::Plane_3` object.
 	*/
	template <typename Traits>
	class Plane_map
	{
	public:
		typedef CGAL::Shape_detection_3::Plane<Traits> Plane_shape;
		typedef boost::shared_ptr<Plane_shape> key_type;
		typedef typename Traits::Plane_3 value_type;
		typedef value_type reference;
		typedef boost::read_write_property_map_tag category;

		inline friend reference get (const Plane_map&, const key_type& k)
		{
			return value_type(*k);
		}

		inline friend void put (const Plane_map&, const key_type& k, const value_type& v)
		{
			k->update(v);
		}
	};
}

}

}

#endif // CGAL_SHAPE_DETECTION_RANSAC_ON_POINTS_3_PLANE_MAP_H
