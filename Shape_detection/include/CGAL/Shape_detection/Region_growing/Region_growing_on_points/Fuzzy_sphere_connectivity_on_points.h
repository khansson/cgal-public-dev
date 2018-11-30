#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_FUZZY_SPHERE_CONNECTIVITY_ON_POINTS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_FUZZY_SPHERE_CONNECTIVITY_ON_POINTS_H

// STL includes.
#include <typeinfo>
#include <type_traits>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>

// Todo:
// Add default values to the constructor.
// Add clear function.

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Fuzzy sphere neighbors search on a set of `Point_2` or `Point_3`.
            \tparam RegionGrowingTraits `CGAL::Shape_detection::Region_growing_traits`
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class RegionGrowingTraits>
        class Fuzzy_sphere_connectivity_on_points {

        public:
        
            using Kernel                  = typename RegionGrowingTraits::Kernel;
            
            using Input_range             = typename RegionGrowingTraits::Input_range;

            using Point                   = typename RegionGrowingTraits::Element;
            ///< Point type, can only be `Point_2` or `Point_3`.

            using Element_map             = typename RegionGrowingTraits::Element_map;
            ///< An `LvaluePropertyMap` that maps to Point.

            using Element_with_properties = typename Element_map::key_type;
            ///< The value type of iterators in `Input_range`.

            #ifdef DOXYGEN_RUNNING
                using Search_base       = unspecified_type;
                ///< Can be `CGAL::Search_traits_2` or `CGAL::Search_traits_3`, automatically deduced based on whether the point type is Point_2 or Point_3.

                using Search_structures = unspecified_type;
                ///< Kd tree configuration class, automatically deduced based on whether `Element_map` is a `CGAL::Identity_property_map` or not. If the `Element_map` is an identity map, it will use the traits class Search_base directly to construct the tree, otherwise the program will create a traits adapter and use that instead.
            #else
                using Search_base       = typename std::conditional<std::is_same<typename Kernel::Point_2, Point>::value, CGAL::Search_traits_2<Kernel>, CGAL::Search_traits_3<Kernel> >::type;
                
				// Here, we have two different specifications to avoid the bug with fuzzy sphere usage in the case when ElementWithProperties = PointType!
				// That is when we have identity property map.
				// Primary template.
                template<class PointType, class ElementWithProperties>
                struct Search_structures {
                    
					using Search_traits_adapter = CGAL::Search_traits_adapter<Element_with_properties, Element_map, Search_base>;
                    using Fuzzy_sphere          = CGAL::Fuzzy_sphere<Search_traits_adapter>;
                    using Tree                  = CGAL::Kd_tree<Search_traits_adapter>;
                };

                // Partial specialization when PointType and ElementWithProperties are the same.
                template<class PointType>
                struct Search_structures<PointType, PointType> {
                    
					using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_base>;
                    using Tree         = CGAL::Kd_tree<Search_base>;
                };
            #endif

            // Element = Point_2 or Point_3.

            using FT             = typename Kernel::FT; 
			///< Number type.

            using Kd_tree_config = Search_structures<Point, Element_with_properties>;
            ///< Kd tree configuration.
                
            using Fuzzy_sphere   = typename Kd_tree_config::Fuzzy_sphere;
            ///< A `CGAL::Fuzzy_sphere` representing the search space of the algorithm.
                
            using Tree           = typename Kd_tree_config::Tree;
            ///< A `CGAL::Kd_tree` holding the points given in the input range.

            /*!
                The constructor takes the point set given in `input_range` and the searching radius, then initializes a Kd tree upon the point set.
            */
            Fuzzy_sphere_connectivity_on_points(const Input_range &input_range, const FT radius) :
            m_input_range(input_range),
            m_radius(radius) {
                
				m_tree.insert(m_input_range.begin(), m_input_range.end());
            }

            /*!
                From a query point `center`, this function creates a `CGAL::Fuzzy_sphere` using the radius previously given in the constructor. It then uses CGAL::Kd_tree::search() to look for the neighbors and push them to `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Element_with_properties &center, Neighbors &neighbors) {

                neighbors.clear();
                Fuzzy_sphere sphere(center, m_radius);
                m_tree.search(std::back_inserter(neighbors), sphere);
            }

        private:

            // Fields.
            const Input_range       &m_input_range;

            const FT                m_radius;
            Tree                    m_tree;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_FUZZY_SPHERE_CONNECTIVITY_ON_POINTS_H
