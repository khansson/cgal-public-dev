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
#include <CGAL/Shape_detection/Region_growing/Tools/Hashable_item_to_index_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Find all faces that share an edge with a given face.
            \tparam FaceGraph General face graph. Model of `FaceGraph`.
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class FaceGraph, 
        class FaceRange = typename FaceGraph::Face_range, 
        class IndexToFaceMap = CGAL::Shape_detection::Random_access_index_to_item_property_map<FaceRange>,
        class FaceToIndexMap = CGAL::Shape_detection::Hashable_item_to_index_property_map<FaceRange> >
        class Connectivity_on_face_graph {

        public:

            using Face_graph                   = FaceGraph;
            ///< General face graph. Model of `FaceGraph`.

            using Face_range                   = FaceRange;
            ///< An `LvaluePropertyMap` that maps to `Face_descriptor`.

            using Face_to_index_map            = FaceToIndexMap;

            using Index_to_face_map            = IndexToFaceMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

            ///< \cond SKIP_IN_MANUAL
            using Index                        = long;
            ///< \endcond

            /*!
                The constructor takes an arbitrary face graph.
                In addition, you can provide an instance of the face descriptor map class.
            */
            Connectivity_on_face_graph(const Face_graph &face_graph) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(m_face_range),
            m_face_to_index_map(m_face_range)
            { }

            Connectivity_on_face_graph(const Face_graph &face_graph, const Index_to_face_map index_to_face_map) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(index_to_face_map),
            m_face_to_index_map(m_face_range)
            { }

            Connectivity_on_face_graph(const Face_graph &face_graph, const Index_to_face_map index_to_face_map, const Face_to_index_map face_to_index_map) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(index_to_face_map),
            m_face_to_index_map(face_to_index_map)
            { }

            /*!
                Using a query index `query_index`, this function retrieves indices of all neighboring faces and stores them in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Index query_index, Neighbors &neighbors) {
                    
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
