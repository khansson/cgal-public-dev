#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_DEFAULT_ITEM_INDEX_TO_ITEM_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_DEFAULT_ITEM_INDEX_TO_ITEM_PROPERTY_MAP_H

namespace CGAL {

    namespace Shape_detection {

        template<class InputRange>
        class Default_item_index_to_item_property_map {
                        
        public:
                    
            using Input_range = InputRange;

            using value_type  = typename InputRange::const_iterator;
            using reference   = const value_type&;
            using key_type    = std::size_t;
            using category    = boost::lvalue_property_map_tag;

            Default_item_index_to_item_property_map(const Input_range &input_range) : 
            m_input_range(input_range)
            { }

            reference operator[](key_type k) const { 
                return m_input_range.begin() + k;
            }

            friend inline reference get(const Default_item_index_to_item_property_map &default_item_index_to_item_map, key_type k) { 
                return default_item_index_to_item_map[k];
            }
                
        private:
            const Input_range &m_input_range;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_DEFAULT_ITEM_INDEX_TO_ITEM_PROPERTY_MAP_H
