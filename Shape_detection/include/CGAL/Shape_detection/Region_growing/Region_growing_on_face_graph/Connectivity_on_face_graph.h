#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_FACE_GRAPH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_FACE_GRAPH_H

// STL includes.
// #include <type_traits>

// CGAL includes.
#include <CGAL/Iterator_range.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Tools/Item_property_map.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Find all faces that share an edge with a given face.
            \tparam FaceGraph General face graph. Model of `FaceGraph`.
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class FaceGraph, class FaceDescriptorMap,
        class IndexToItemMap = CGAL::Shape_detection::Random_access_index_to_item_property_map<typename FaceGraph::Face_range> >
        class Connectivity_on_face_graph {

        public:

            using Face_graph                   = FaceGraph;
            ///< General face graph. Model of `FaceGraph`.

            using Face_descriptor_map          = FaceDescriptorMap;
            ///< An `LvaluePropertyMap` that maps to `Face_descriptor`.

            using Index_to_item_map            = IndexToItemMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

            ///< \cond SKIP_IN_MANUAL
            using Index                        = std::size_t;

            using Index_to_face_descriptor_map = CGAL::Shape_detection::Item_property_map<Index_to_item_map, Face_descriptor_map>;
            ///< \endcond

            /*!
                The constructor takes an arbitrary face graph.
                In addition, you can provide an instance of the face descriptor map class.
            */
            Connectivity_on_face_graph(const Face_graph &face_graph, const Face_descriptor_map &face_descriptor_map = Face_descriptor_map()) :
            m_face_graph(face_graph),
            m_index_to_item_map(m_face_graph.faces()),
            m_index_to_face_descriptor_map(m_index_to_item_map, face_descriptor_map) 
            { }

            Connectivity_on_face_graph(const Face_graph &face_graph, const Index_to_item_map &index_to_item_map, const Face_descriptor_map &face_descriptor_map = Face_descriptor_map()) :
            m_face_graph(face_graph),
            m_index_to_item_map(index_to_item_map),
            m_index_to_face_descriptor_map(m_index_to_item_map, face_descriptor_map) 
            { }

            /*!
                Using a query index `query_index`, this function retrieves indices of all neighboring faces and stores them in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Index query_index, Neighbors &neighbors) {
                    
                neighbors.clear();

                const auto &query_face_descriptor     = get(m_index_to_face_descriptor_map, query_index);
                const auto &query_halfedge_descriptor = CGAL::halfedge(query_face_descriptor, m_face_graph);

                const auto &faces = CGAL::faces_around_face(query_halfedge_descriptor, m_face_graph);

                for (auto face = faces.begin(); face != faces.end(); ++face) {
                    if (*face != m_face_graph.null_face()) {
                        
                        const Index face_index = static_cast<Index>(*face);
                        neighbors.push_back(face_index);
                    }
                }
            }

        private:

            // Fields.
            const Face_graph                       &m_face_graph;

            const Index_to_item_map                 m_index_to_item_map;
            const Index_to_face_descriptor_map      m_index_to_face_descriptor_map;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_FACE_GRAPH_H
