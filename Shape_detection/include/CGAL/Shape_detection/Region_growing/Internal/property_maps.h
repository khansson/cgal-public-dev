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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_PROPERTY_MAPS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_PROPERTY_MAPS_H

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<class ItemRange, class PropertyMap>
  class Item_property_map {
                        
  public: 
    using Item_range = ItemRange;
    using Property_map = PropertyMap;

    using value_type = typename Property_map::value_type;
    using reference = const value_type&;
    using key_type = std::size_t;
    using category = boost::lvalue_property_map_tag;

    Item_property_map(
      const Item_range& item_range, 
      const Property_map& property_map) : 
    m_item_range(item_range),
    m_property_map(property_map) 
    { }

    reference operator[](key_type item_index) const { 
                
      const auto& key = *(m_item_range.begin() + item_index);
      return get(m_property_map, key);
    }

    friend inline reference get(
      const Item_property_map& item_map, 
      key_type key) { 
      
      return item_map[key];
    }
                
  private:
    const Item_range& m_item_range;
    const Property_map& m_property_map;
  };

} // internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_PROPERTY_MAPS_H
