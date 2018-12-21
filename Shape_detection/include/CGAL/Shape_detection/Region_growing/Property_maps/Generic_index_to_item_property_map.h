#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_GENERIC_INDEX_TO_ITEM_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_GENERIC_INDEX_TO_ITEM_PROPERTY_MAP_H

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {

    namespace Shape_detection {

        template<class InputRange>
        class Generic_index_to_item_property_map {
                        
        public:
                    
            using Input_range = InputRange;

            using value_type  = typename InputRange::const_iterator;
            using key_type    = long;
            using category    = boost::lvalue_property_map_tag;

            Generic_index_to_item_property_map(const Input_range &input_range) : 
            m_input_range(input_range)
            { }

            value_type operator[](key_type k) const { 
                
                key_type count = 0;
                CGAL_precondition(k < m_input_range.size());
                
                for (value_type item = m_input_range.begin(); item != m_input_range.end(); ++item, ++count)
                    if (count == k) return item;
            }

            friend inline value_type get(const Generic_index_to_item_property_map &index_to_item_map, key_type k) { 
                return index_to_item_map[k];
            }
                
        private:
            const Input_range &m_input_range;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_GENERIC_INDEX_TO_ITEM_PROPERTY_MAP_H
