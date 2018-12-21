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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_FACE_GRAPH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_FACE_GRAPH_H

// Boost includes.
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Face graph includes.
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// CGAL includes.
#include <CGAL/assertions.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Property_maps/Hashable_item_to_index_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Property_maps/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Find all faces that share an edge with a given face.
            \tparam FaceGraph General face graph. Model of `FaceGraph`.
            \tparam FaceRange An arbitrary range with graph faces, given an IndexToFaceMap is provided. The default one is random access.
            \tparam IndexToFaceMap An `LvaluePropertyMap` that maps face to `Index`, which is any signed integer type, the default `Index` is `long`.
            \tparam FaceToIndexMap An `LvaluePropertyMap` that `Index` to face, opposite to IndexToFaceMap.
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class FaceGraph, 
        class FaceRange = typename FaceGraph::Face_range, 
        class IndexToFaceMap = CGAL::Shape_detection::Random_access_index_to_item_property_map<FaceRange>,
        class FaceToIndexMap = CGAL::Shape_detection::Hashable_item_to_index_property_map<FaceRange> >
        class Connectivity_on_face_graph {

        public:

            /// \name Types
            /// @{

            using Face_graph                   = FaceGraph;
            ///< General face graph. Model of `FaceGraph`.

            using Face_range                   = FaceRange;
            ///< An arbitrary range with graph faces.

            using Index_to_face_map            = IndexToFaceMap;
            ///< An `LvaluePropertyMap` that maps `Index` to face.

            using Face_to_index_map            = FaceToIndexMap;
            ///< An `LvaluePropertyMap` that maps face index to face.

            ///< \cond SKIP_IN_MANUAL
            using Index                        = long;
            ///< \endcond

            /// @}

            /// \name Initialization
            /// @{

            /*!
                The constructor that takes an arbitrary face graph.
            */
            Connectivity_on_face_graph(const Face_graph &face_graph) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(m_face_range),
            m_face_to_index_map(m_face_range)
            { }

            /*!
                The constructor that takes an arbitrary face graph and an index_to_face_map to access a face given its index in the face range from the `face_graph`.
            */
            Connectivity_on_face_graph(const Face_graph &face_graph, const Index_to_face_map index_to_face_map) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(index_to_face_map),
            m_face_to_index_map(m_face_range)
            { }

            /*!
                The constructor that takes an arbitrary face graph, an index_to_face_map to access a face given its index in the face range from the `face_graph`,
                and a face_to_index_map to access face index given a face.
            */
            Connectivity_on_face_graph(const Face_graph &face_graph, const Index_to_face_map index_to_face_map, const Face_to_index_map face_to_index_map) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(index_to_face_map),
            m_face_to_index_map(face_to_index_map)
            { }

            /// @}

            /// \name Access
            /// @{ 

            /*!
                Using a query index `query_index`, this function retrieves indices of all neighboring faces and stores them in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Index query_index, Neighbors &neighbors) const {
                    
                CGAL_precondition(query_index < m_face_range.size());
                neighbors.clear();

                const auto &query_face     = get(m_index_to_face_map, query_index);
                const auto &query_halfedge = CGAL::halfedge(*query_face, m_face_graph);

                const auto &faces = CGAL::faces_around_face(query_halfedge, m_face_graph);
                for (auto face = faces.begin(); face != faces.end(); ++face) {

                    const Index face_index = get(m_face_to_index_map, *face);
                    if (face_index >= 0) neighbors.push_back(face_index);
                }
            }

            /// @}

        private:

            // Fields.
            const Face_graph           &m_face_graph;
            const Face_range            m_face_range;

            const Index_to_face_map     m_index_to_face_map;
            const Face_to_index_map     m_face_to_index_map;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_FACE_GRAPH_H
