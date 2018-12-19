#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_SQRT_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_SQRT_H

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/number_utils.h>

namespace CGAL {

    namespace Shape_detection {

        template<class Traits> 
        class Default_sqrt {
            using FT = typename Traits::FT;

        public:
            FT operator()(const FT value) const { 
                return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
            }
        };

        BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

        // Case: do_not_use_default = false.
        template<class Traits, bool do_not_use_default = Has_nested_type_Sqrt<Traits>::value>
        class Get_sqrt {
        
        public:
            using Sqrt = Default_sqrt<Traits>;

            static Sqrt sqrt_object(const Traits &) { 
                return Sqrt();
            }
        };

        // Case: do_not_use_default = true.
        template<class Traits>
        class Get_sqrt<Traits, true> {
        
        public:
            using Sqrt = typename Traits::Sqrt;

            static Sqrt sqrt_object(const Traits &traits) { 
                return traits.sqrt_object();
            }
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_SQRT_H
