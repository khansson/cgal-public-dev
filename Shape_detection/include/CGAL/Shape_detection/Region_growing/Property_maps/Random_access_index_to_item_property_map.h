// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_RANDOM_ACCESS_INDEX_TO_ITEM_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_RANDOM_ACCESS_INDEX_TO_ITEM_PROPERTY_MAP_H

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {

    namespace Shape_detection {

        /*! 
            \ingroup PkgShapeDetection
            \brief An `LvaluePropertyMap` that maps `Index` to an item from the InputRange, where `Index` is any signed integer type, the default `Index` is `long`. This map works only on a random access container.
            \tparam InputRange A random access range with user-defined items.
        */
        template<class InputRange>
        class Random_access_index_to_item_property_map {
                        
        public:
            
            /// \name Types
            /// @{

            using Input_range = InputRange;
            ///< A random access range with user-defined items.

            ///< \cond SKIP_IN_MANUAL
            using value_type  = typename InputRange::const_iterator;
            using key_type    = long;
            using category    = boost::lvalue_property_map_tag;
            ///< \endcond

            /// @}

            /// \name Initialization
            /// @{

            /*!
                The constructor that takes a set of items stored in a random access container `input_range`.
            */
            Random_access_index_to_item_property_map(const Input_range &input_range) : 
            m_input_range(input_range)
            { }

            /// @}

            ///< \cond SKIP_IN_MANUAL
            value_type operator[](key_type k) const { 

                CGAL_precondition(k < m_input_range.size());
                return m_input_range.begin() + k;
            }

            friend inline value_type get(const Random_access_index_to_item_property_map &index_to_item_map, key_type k) { 
                return index_to_item_map[k];
            }
            ///< \endcond
                
        private:
            const Input_range &m_input_range;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_RANDOM_ACCESS_INDEX_TO_ITEM_PROPERTY_MAP_H
