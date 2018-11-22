#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_TRAITS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_TRAITS_H

namespace CGAL {

    namespace Shape_detection {

            /*!
                \ingroup PkgShapeDetection
                \brief Traits class describing the input values and the primary type of geometric element on which we grow regions.
                \tparam InputRange
                \tparam ElementMap
                \tparam GeomTraits
            */

            template<class InputRange, class ElementMap, class GeomTraits>
            class Region_growing_traits {
            public:

                using Input_range             = InputRange;
                ///< An `Iterator_range` of bidirectional iterator.

                using Element_map             = ElementMap;
                ///< A model of `LvaluePropertyMap` that maps to a geometric type.

                using Element                 = typename Element_map::value_type;
                ///< The type of an element in the given `Input_range` of elements to grow regions.
                
                // Element_map::value_type must be Kernel::Point_2 or Kernel::Point_3 or Face_index.
                // Element_map::key_type must be the same as std::iterator_traits<Input_range::iterator>::value_type.

                using Kernel                  = GeomTraits;
                ///< The kernel with basic geometric types and operations.
            };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_TRAITS_H
