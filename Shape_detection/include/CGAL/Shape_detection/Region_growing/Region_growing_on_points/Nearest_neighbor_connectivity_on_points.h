// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

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
#include <CGAL/Shape_detection/Region_growing/Internal/Item_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Property_maps/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetectionRGOnPoints
            \brief K nearest neighbors (kNN) search on a set of `Point_2` or `Point_3`.
            \tparam Traits Model of `Kernel`
            \tparam InputRange An arbitrary range with user-defined items, given an IndexToItemMap is provided. The default one is random access.
            \tparam PointMap An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.
            \tparam IndexToItemMap An `LvaluePropertyMap` that maps item to `Item_index`, which is any signed integer type, the default `Item_index` is `long`.
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class Traits, class InputRange, class PointMap,
        class IndexToItemMap = CGAL::Shape_detection::Random_access_index_to_item_property_map<InputRange> >
        class Nearest_neighbor_connectivity_on_points {

        public:

            /// \name Types
            /// @{
            
            using Input_range             = InputRange;
            ///< An arbitrary range with user-defined items.

            using Point_map               = PointMap;
            ///< An `LvaluePropertyMap` that maps to `Point_2` or `Point_3`.

            using Index_to_item_map       = IndexToItemMap;
            ///< An `LvaluePropertyMap` that maps `Item_index` to item.

            using Point                   = typename Point_map::value_type;
            ///< Point type, can only be `Point_2` or `Point_3`.

            ///< \cond SKIP_IN_MANUAL
            using Item_index              = long;

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
                    
					using Search_traits   = CGAL::Search_traits_adapter<Item_index, Index_to_point_map, Search_base>;
                    using Distance        = CGAL::Distance_adapter<Item_index, Index_to_point_map, CGAL::Euclidean_distance<Search_base> >;
                    using Splitter        = CGAL::Sliding_midpoint<Search_traits>;
                    using Search_tree     = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;
                    using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Search_tree>;
                    using Tree            = typename Neighbor_search::Tree;
                };

            #endif
                
            using Neighbor_search         = typename Search_structures::Neighbor_search;
            ///< A search class `CGAL::Orthogonal_k_neighbor_search` that implements a Kd tree.
                
            ///< \cond SKIP_IN_MANUAL
            using Distance                = typename Search_structures::Distance;
            ///< \endcond

            using Tree                    = typename Search_structures::Tree;
            ///< The Kd tree member type of the search class.

            /// @}

            /// \name Initialization
            /// @{

            /*!
                The constructor that takes a set of items, provided a point_map to access a `Point` from an item, and a number of nearest neighbors (the value "k" in "kNN").
            */
            Nearest_neighbor_connectivity_on_points(const Input_range &input_range, const std::size_t number_of_neighbors = 12, const Point_map point_map = Point_map()) :
            m_input_range(input_range),
            m_number_of_neighbors(number_of_neighbors),
            m_index_to_item_map(m_input_range),
            m_point_map(point_map),
            m_index_to_point_map(m_index_to_item_map, m_point_map),
            m_distance(m_index_to_point_map),
            m_tree(
                boost::counting_iterator<Item_index>(0),
                boost::counting_iterator<Item_index>(m_input_range.size()),
                typename Search_structures::Splitter(),
                typename Search_structures::Search_traits(m_index_to_point_map)) { 

                    m_tree.build();
                    CGAL_precondition(number_of_neighbors >= 0);
                }

            /*!
                The constructor that takes a set of items, provided a point_map to access a `Point` from an item, a number of nearest neighbors (the value "k" in "kNN"), and an index_to_item_map to access an item given its index in the `input_range`.
            */
            Nearest_neighbor_connectivity_on_points(const Input_range &input_range, const Index_to_item_map index_to_item_map, const std::size_t number_of_neighbors = 12, const Point_map point_map = Point_map()) :
            m_input_range(input_range),
            m_number_of_neighbors(number_of_neighbors),
            m_index_to_item_map(index_to_item_map),
            m_point_map(point_map),
            m_index_to_point_map(m_index_to_item_map, m_point_map),
            m_distance(m_index_to_point_map),
            m_tree(
                boost::counting_iterator<Item_index>(0),
                boost::counting_iterator<Item_index>(m_input_range.size()),
                typename Search_structures::Splitter(),
                typename Search_structures::Search_traits(m_index_to_point_map)) { 

                    m_tree.build();
                    CGAL_precondition(number_of_neighbors >= 0);
                }

            /// @}

            /// \name Access
            /// @{ 

            /*!
                This function takes the index `query_index` of a query item and returns indices of the k closest items around it. The result is stored in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Item_index query_index, Neighbors &neighbors) const {
                CGAL_precondition(query_index < m_input_range.size());

                neighbors.clear();
                Neighbor_search neighbor_search(m_tree, get(m_index_to_point_map, query_index), m_number_of_neighbors, 0, true, m_distance);
                
                for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
                    neighbors.push_back(it->first);
            }

            /// @}

            /// \name Memory Management
            /// @{

            /*!
                Clear all internal data structures.
            */
            void clear() {
                m_tree.clear();
            }

            /// @}

        private:

            // Fields.
            const Input_range              &m_input_range;
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
