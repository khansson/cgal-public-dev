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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_HASHABLE_ITEM_TO_INDEX_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_HASHABLE_ITEM_TO_INDEX_PROPERTY_MAP_H

// STL includes.
#include <unordered_map>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {

    namespace Shape_detection {

        /*! 
            \ingroup PkgShapeDetection
            \brief An `LvaluePropertyMap` that maps a hashable item to `Index`, where `Index` is any signed integer type, the default `Index` is `long`. This map works on an arbitrary container with items.
            \tparam InputRange A range with user-defined items, where each item is hashable.
        */
        template<class InputRange>
        class Hashable_item_to_index_property_map {
                        
        public:
                    
            /// \name Types
            /// @{

            using Input_range = InputRange;
            ///< A range with user-defined items. Each item is hashable and stored in an arbitrary container.
            
            ///< \cond SKIP_IN_MANUAL
            using Iterator = typename Input_range::const_iterator;
            using Item     = typename Iterator::value_type;

            using value_type = long;
            using key_type   = Item;
            using category   = boost::lvalue_property_map_tag;

            using Item_map = std::unordered_map<key_type, value_type>;
            ///< \endcond

            /// @}

            /// \name Initialization
            /// @{

            /*!
                The constructor that takes a set of hashable items stored in an arbitrary container `input_range`.
            */
            Hashable_item_to_index_property_map(const Input_range &input_range) : 
            m_input_range(input_range) { 

                value_type i = 0;
                for (Iterator item = m_input_range.begin(); item != m_input_range.end(); ++item, ++i)
                    m_item_map[*item] = i;
            }

            /// @}

            ///< \cond SKIP_IN_MANUAL
            value_type operator[](const key_type &k) const { 

                const auto &value = m_item_map.find(k);
                if (value == m_item_map.end()) 
                    return -1;
                return value->second;
            }

            friend inline value_type get(const Hashable_item_to_index_property_map &item_to_index_map, const key_type &k) { 
                return item_to_index_map[k];
            }
            ///< \endcond
                
        private:
            const Input_range &m_input_range;
            Item_map m_item_map;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_HASHABLE_ITEM_TO_INDEX_PROPERTY_MAP_H
