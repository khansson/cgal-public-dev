#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_NEAREST_NEIGHBOR_CONNECTIVITY_ON_POINTS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_NEAREST_NEIGHBOR_CONNECTIVITY_ON_POINTS_H

// STL includes.
#include <typeinfo>
#include <type_traits>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief K nearest neighbors (kNN) search on a set of `Point_2` or `Point_3`.
            \tparam RegionGrowingTraits `CGAL::Shape_detection::Region_growing_traits`
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class RegionGrowingTraits>
        class Nearest_neighbor_connectivity_on_points {

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

                // Here, we have two different specifications to avoid the bug with k nearest neighbor usage in the case when ElementWithProperties = PointType!
				// That is when we have identity property map.
				// Primary template.
                template<class PointType, class ElementWithProperties>
                struct Search_structures {
                    
                    using Search_traits_adapter = CGAL::Search_traits_adapter<Element_with_properties, Element_map, Search_base>;
                    using Neighbor_search       = CGAL::Orthogonal_k_neighbor_search<Search_traits_adapter>;
                    using Tree                  = typename Neighbor_search::Tree;
                };

                // Partial specialization when PointType and ElementWithProperties are the same.
                template<class PointType>
                struct Search_structures<PointType, PointType> {

                    using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_base>;
                    using Tree            = typename Neighbor_search::Tree;
                };
            #endif

            // Element = Point_2 or Point_3.

            using FT              = typename Kernel::FT;
            ///< Number type.

            using Kd_tree_config  = Search_structures<Point, Element_with_properties>;
            ///< Kd tree configuration.

            using Neighbor_search = typename Kd_tree_config::Neighbor_search;
            ///< A search class CGAL::Orthogonal_k_neighbor_search that implements a Kd tree.

            using Tree            = typename Kd_tree_config::Tree;
            ///< The Kd tree member type of the search class.

            /*!
                The constructor initializes a Kd tree using the user-provided input and stores `number_of_neighbors` (the value of "k", as in "kNN") for later use.
            */
            Nearest_neighbor_connectivity_on_points(const Input_range &input_range, const size_t number_of_neighbors) :
            m_input_range(input_range),
            m_number_of_neighbors(number_of_neighbors) {
            
                m_tree.insert(m_input_range.begin(), m_input_range.end());
            }

            /*!
                The function takes a query point and return the K closest points around it. The result is stored in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Element_with_properties &query, Neighbors &neighbors) {

                neighbors.clear();
                Neighbor_search search(m_tree, get(m_elem_map, query), m_number_of_neighbors);
                for (auto it = search.begin(); it != search.end(); ++it)
                    neighbors.push_back(it->first);
            }

        private:

            // Fields.
            const Input_range       &m_input_range;

            const size_t             m_number_of_neighbors;
            const Element_map        m_elem_map = Element_map();
            Tree                     m_tree;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_NEAREST_NEIGHBOR_CONNECTIVITY_ON_POINTS_H
