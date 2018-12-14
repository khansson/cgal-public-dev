#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_FUZZY_SPHERE_CONNECTIVITY_ON_POINTS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_FUZZY_SPHERE_CONNECTIVITY_ON_POINTS_H

// STL includes.
#include <vector>
#include <typeinfo>
#include <type_traits>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/assertions.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Tools/Item_to_point_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Default_item_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Fuzzy sphere search for neighbors on a set of `Point_2` or `Point_3`.
            \tparam InputRange A random access range with user-defined items.
            \tparam PointMap An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.
            \tparam Traits Model of `RegionGrowingOnPointsTraits`
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class InputRange, class PointMap, class Traits,
        class ItemIndexToItemMap = CGAL::Shape_detection::Default_item_index_to_item_property_map<InputRange> >
        class Fuzzy_sphere_connectivity_on_points {

        public:
            
            using Input_range             = InputRange;
            ///< An arbitrary range with user-defined items.

            using Point_map               = PointMap;
            ///< An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.

            using Point                   = typename PointMap::value_type;
            ///< Point type, can only be `Point_2` or `Point_3`.

            using Item_index_to_item_map  = ItemIndexToItemMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

            ///< \cond SKIP_IN_MANUAL
            using Item_index              = std::size_t;

            using Item_index_to_point_map = CGAL::Shape_detection::Item_to_point_property_map<Item_index_to_item_map, Point_map>;
            ///< \endcond

            #ifdef DOXYGEN_RUNNING
                
                using Search_base         = unspecified_type;
                ///< Can be `CGAL::Search_traits_2` or `CGAL::Search_traits_3`, automatically deduced based on whether the point type is `Point_2` or `Point_3`.

                using Search_structures   = unspecified_type;
                ///< Kd tree configuration class.

            #else
                
                using Search_base         = typename std::conditional<std::is_same<typename Traits::Point_2, Point>::value, CGAL::Search_traits_2<Traits>, CGAL::Search_traits_3<Traits> >::type;

                struct Search_structures {
                    
					using Search_traits   = CGAL::Search_traits_adapter<Item_index, Item_index_to_point_map, Search_base>;
                    using Splitter        = CGAL::Sliding_midpoint<Search_traits>;
                    using Fuzzy_sphere    = CGAL::Fuzzy_sphere<Search_traits>;
                    using Tree            = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;
                };

            #endif

            using FT                      = typename Traits::FT;
			///< Number type.
                
            using Fuzzy_sphere            = typename Search_structures::Fuzzy_sphere;
            ///< Represents the search space of the algorithm. It is defined as `CGAL::Fuzzy_sphere`.
                
            using Tree                    = typename Search_structures::Tree;
            ///< `CGAL::Kd_tree` over the items from the input range.

            /*!
                The constructor takes the point set given in `input_range` and the searching radius, then initializes a Kd tree upon the point set.
                In addition, you can provide an instance of the point map class.
            */
            Fuzzy_sphere_connectivity_on_points(const Input_range &input_range, const FT search_radius = FT(1), const Point_map &point_map = PointMap()) :
            m_search_radius(search_radius),
            m_item_index_to_item_map(input_range),
            m_item_index_to_point_map(m_item_index_to_item_map, point_map),
            m_tree(
                boost::counting_iterator<Item_index>(0),
                boost::counting_iterator<Item_index>(input_range.size()),
                typename Search_structures::Splitter(),
                typename Search_structures::Search_traits(m_item_index_to_point_map)) { 

                    m_tree.build();
                    CGAL_precondition(search_radius >= FT(0));
                }

            Fuzzy_sphere_connectivity_on_points(const Input_range &input_range, const Item_index_to_item_map &item_index_to_item_map, const FT search_radius = FT(1), const Point_map &point_map = PointMap()) :
            m_search_radius(search_radius),
            m_item_index_to_item_map(item_index_to_item_map),
            m_item_index_to_point_map(m_item_index_to_item_map, point_map),
            m_tree(
                boost::counting_iterator<Item_index>(0),
                boost::counting_iterator<Item_index>(input_range.size()),
                typename Search_structures::Splitter(),
                typename Search_structures::Search_traits(m_item_index_to_point_map)) { 

                    m_tree.build();
                    CGAL_precondition(search_radius >= FT(0));
                }

            /*!
                From a query item with the index `query_index`, this function creates a search sphere centered at this item.
                It then uses CGAL::Kd_tree::search() to look for the neighbors of the given query and push their indices to `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Item_index query_index, Neighbors &neighbors) const {
                
                neighbors.clear();
                const Fuzzy_sphere sphere(query_index, m_search_radius, FT(0), m_tree.traits());
                m_tree.search(std::back_inserter(neighbors), sphere);
            }

            void clear() {
                m_tree.clear();
            }

        private:

            // Fields.
            const FT                            m_search_radius;
            const Item_index_to_item_map        m_item_index_to_item_map;
            const Item_index_to_point_map       m_item_index_to_point_map;
            Tree                                m_tree;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_FUZZY_SPHERE_CONNECTIVITY_ON_POINTS_H
