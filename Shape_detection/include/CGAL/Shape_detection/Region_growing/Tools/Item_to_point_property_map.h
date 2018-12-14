#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_ITEM_TO_POINT_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_ITEM_TO_POINT_PROPERTY_MAP_H

namespace CGAL {

    namespace Shape_detection {

        template<class ItemIndexToItemMap, class PointMap>
        class Item_to_point_property_map {
                        
        public:
                    
            using Item_index_to_item_map = ItemIndexToItemMap;
            using Point_map              = PointMap;

            using value_type = typename Point_map::value_type;
            using reference  = const value_type&;
            using key_type   = std::size_t;
            using category   = boost::lvalue_property_map_tag;

            Item_to_point_property_map(const Item_index_to_item_map &item_index_to_item_map, const Point_map &point_map) : 
            m_item_index_to_item_map(item_index_to_item_map),
            m_point_map(point_map) 
            { }

            reference operator[](key_type k) const { 
                
                const auto &item = get(m_item_index_to_item_map, k);
                return get(m_point_map, *item);
            }

            friend inline reference get(const Item_to_point_property_map &item_to_point_map, key_type k) { 
                return item_to_point_map[k];
            }
                
        private:
            const Item_index_to_item_map &m_item_index_to_item_map;
            const Point_map              &m_point_map;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_ITEM_TO_POINT_PROPERTY_MAP_H
