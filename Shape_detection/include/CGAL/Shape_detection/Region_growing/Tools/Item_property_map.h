#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_ITEM_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_ITEM_PROPERTY_MAP_H

namespace CGAL {

    namespace Shape_detection {

        template<class IndexToItemMap, class PropertyMap>
        class Item_property_map {
                        
        public:
                    
            using Index_to_item_map = IndexToItemMap;
            using Property_map      = PropertyMap;

            using value_type = typename Property_map::value_type;
            using reference  = const value_type&;
            using key_type   = std::size_t;
            using category   = boost::lvalue_property_map_tag;

            Item_property_map(const Index_to_item_map &index_to_item_map, const Property_map &property_map) : 
            m_index_to_item_map(index_to_item_map),
            m_property_map(property_map) 
            { }

            reference operator[](key_type index) const { 
                
                const auto &item = get(m_index_to_item_map, index);
                return get(m_property_map, *item);
            }

            friend inline reference get(const Item_property_map &item_map, key_type k) { 
                return item_map[k];
            }
                
        private:
            const Index_to_item_map &m_index_to_item_map;
            const Property_map      &m_property_map;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_ITEM_PROPERTY_MAP_H
