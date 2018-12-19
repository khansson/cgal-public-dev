#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_HASHABLE_ITEM_TO_INDEX_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_HASHABLE_ITEM_TO_INDEX_PROPERTY_MAP_H

// STL includes.
#include <unordered_map>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {

    namespace Shape_detection {

        template<class InputRange>
        class Hashable_item_to_index_property_map {
                        
        public:
                    
            using Input_range = InputRange;
            
            using Iterator = typename Input_range::const_iterator;
            using Item     = typename Iterator::value_type;

            using value_type = int;
            using key_type   = Item;
            using category   = boost::lvalue_property_map_tag;

            using Item_map = std::unordered_map<key_type, value_type>;

            Hashable_item_to_index_property_map(const Input_range &input_range) : 
            m_input_range(input_range) { 

                value_type i = 0;
                for (Iterator item = m_input_range.begin(); item != m_input_range.end(); ++item, ++i)
                    m_item_map[*item] = i;
            }

            value_type operator[](const key_type &k) const { 

                const auto &value = m_item_map.find(k);
                if (value == m_item_map.end()) 
                    return -1;
                return value->second;
            }

            friend inline value_type get(const Hashable_item_to_index_property_map &item_to_index_map, const key_type &k) { 
                return item_to_index_map[k];
            }
                
        private:
            const Input_range &m_input_range;
            Item_map m_item_map;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_HASHABLE_ITEM_TO_INDEX_PROPERTY_MAP_H
