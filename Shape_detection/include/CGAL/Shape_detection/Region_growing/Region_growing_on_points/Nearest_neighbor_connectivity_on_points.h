#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_NEAREST_NEIGHBOR_CONNECTIVITY_ON_POINTS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_NEAREST_NEIGHBOR_CONNECTIVITY_ON_POINTS_H

// STL includes.
#include <vector>
#include <typeinfo>
#include <type_traits>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Tools/Item_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief K nearest neighbors (kNN) search on a set of `Point_2` or `Point_3`.
            \tparam InputRange A random access range with user-defined items.
            \tparam PointMap An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.
            \tparam Traits Model of `RegionGrowingOnPointsTraits`
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class InputRange, class PointMap, class Traits,
        class IndexToItemMap = CGAL::Shape_detection::Random_access_index_to_item_property_map<InputRange> >
        class Nearest_neighbor_connectivity_on_points {

        public:
            
            using Input_range             = InputRange;
            ///< An arbitrary range with user-defined items.

            using Point_map               = PointMap;
            ///< An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.

            using Point                   = typename Point_map::value_type;
            ///< Point type, can only be `Point_2` or `Point_3`.

            using Index_to_item_map       = IndexToItemMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

            ///< \cond SKIP_IN_MANUAL
            using Index                   = std::size_t;

            using Index_to_point_map      = CGAL::Shape_detection::Item_property_map<Index_to_item_map, Point_map>;
            ///< \endcond

            #ifdef DOXYGEN_RUNNING
                
                using Search_base         = unspecified_type;
                ///< Can be `CGAL::Search_traits_2` or `CGAL::Search_traits_3`, automatically deduced based on whether the point type is `Point_2` or `Point_3`.

                using Search_structures   = unspecified_type;
                ///< Kd tree configuration class.

            #else
                
                using Search_base         = typename std::conditional<std::is_same<typename Traits::Point_2, Point>::value, CGAL::Search_traits_2<Traits>, CGAL::Search_traits_3<Traits> >::type;

                struct Search_structures {
                    
					using Search_traits   = CGAL::Search_traits_adapter<Index, Index_to_point_map, Search_base>;
                    using Distance        = CGAL::Distance_adapter<Index, Index_to_point_map, CGAL::Euclidean_distance<Search_base> >;
                    using Splitter        = CGAL::Sliding_midpoint<Search_traits>;
                    using Search_tree     = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;
                    using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Search_tree>;
                    using Tree            = typename Neighbor_search::Tree;
                };

            #endif
                
            using Neighbor_search         = typename Search_structures::Neighbor_search;
            ///< A search class CGAL::Orthogonal_k_neighbor_search that implements a Kd tree.
                
            using Distance                = typename Search_structures::Distance;
            ///< A distance class for a Kd tree.

            using Tree                    = typename Search_structures::Tree;
            ///< The Kd tree member type of the search class.

            /*!
                The constructor takes the point set given in `input_range` and the number of nearest neighbors (the value "k" in "kNN"), then initializes a Kd tree upon the point set.
                In addition, you can provide an instance of the point map class.
            */
            Nearest_neighbor_connectivity_on_points(const Input_range &input_range, const std::size_t number_of_neighbors = 12, const Point_map point_map = Point_map()) :
            m_number_of_neighbors(number_of_neighbors),
            m_index_to_item_map(input_range),
            m_point_map(point_map),
            m_index_to_point_map(m_index_to_item_map, m_point_map),
            m_distance(m_index_to_point_map),
            m_tree(
                boost::counting_iterator<Index>(0),
                boost::counting_iterator<Index>(input_range.size()),
                typename Search_structures::Splitter(),
                typename Search_structures::Search_traits(m_index_to_point_map)) { 

                    m_tree.build();
                    CGAL_precondition(number_of_neighbors >= 0);
                }

            Nearest_neighbor_connectivity_on_points(const Input_range &input_range, const Index_to_item_map index_to_item_map, const std::size_t number_of_neighbors = 12, const Point_map point_map = Point_map()) :
            m_number_of_neighbors(number_of_neighbors),
            m_index_to_item_map(index_to_item_map),
            m_point_map(point_map),
            m_index_to_point_map(m_index_to_item_map, m_point_map),
            m_distance(m_index_to_point_map),
            m_tree(
                boost::counting_iterator<Index>(0),
                boost::counting_iterator<Index>(input_range.size()),
                typename Search_structures::Splitter(),
                typename Search_structures::Search_traits(m_index_to_point_map)) { 

                    m_tree.build();
                    CGAL_precondition(number_of_neighbors >= 0);
                }

            /*!
                This function takes an index of a query item and return indices of the k closest items around it. The result is stored in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Index query_index, Neighbors &neighbors) const {
                
                neighbors.clear();
                Neighbor_search neighbor_search(m_tree, get(m_index_to_point_map, query_index), m_number_of_neighbors, 0, true, m_distance);
                
                for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
                    neighbors.push_back(it->first);
            }

            void clear() {
                m_tree.clear();
            }

        private:

            // Fields.
            const std::size_t               m_number_of_neighbors;

            const Index_to_item_map         m_index_to_item_map;
            const Point_map                 m_point_map;
            const Index_to_point_map        m_index_to_point_map;

            Distance                        m_distance;
            Tree                            m_tree;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_NEAREST_NEIGHBOR_CONNECTIVITY_ON_POINTS_H
