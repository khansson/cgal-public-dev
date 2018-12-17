#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_RANDOM_ACCESS_INDEX_TO_ITEM_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_RANDOM_ACCESS_INDEX_TO_ITEM_PROPERTY_MAP_H

namespace CGAL {

    namespace Shape_detection {

        template<class InputRange>
        class Random_access_index_to_item_property_map {
                        
        public:
                    
            using Input_range = InputRange;

            using value_type  = typename InputRange::const_iterator;
            using key_type    = std::size_t;
            using category    = boost::lvalue_property_map_tag;

            Random_access_index_to_item_property_map(const Input_range &input_range) : 
            m_input_range(input_range)
            { }

            value_type operator[](key_type k) const { 
                return static_cast<value_type>(m_input_range.begin() + k);
            }

            friend inline value_type get(const Random_access_index_to_item_property_map &index_to_item_map, key_type k) { 
                return index_to_item_map[k];
            }
                
        private:
            const Input_range &m_input_range;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_RANDOM_ACCESS_INDEX_TO_ITEM_PROPERTY_MAP_H
