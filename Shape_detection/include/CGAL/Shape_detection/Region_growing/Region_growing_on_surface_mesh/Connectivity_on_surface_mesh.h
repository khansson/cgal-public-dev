#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_SURFACE_MESH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_SURFACE_MESH_H

// CGAL includes.
#include <CGAL/Iterator_range.h>

// Todo:
// Add const to Face face.
// Remove Face(4294967295).
// Can we generalize this class wrt face graph concept?

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Find all faces that share an edge with a given face.
            \tparam RegionGrowingTraits `CGAL::Shape_detection::Region_growing_traits`
            \tparam Mesh Has model `CGAL::Surface_mesh`.
            \cgalModels `RegionGrowingConnectivity`
        */
        template<class RegionGrowingTraits, class Mesh>
        class Connectivity_on_surface_mesh {

        public:
        
            using Kernel                  = typename RegionGrowingTraits::Kernel;

            using Face_range              = CGAL::Iterator_range< CGAL::Face_around_face_iterator<Mesh> >;
            ///< Result type of CGAL::faces_around_face().

            using Face                    = typename RegionGrowingTraits::Element;
            ///< Must be equivalent to Mesh::Face_index.
                
            using Element_map             = typename RegionGrowingTraits::Element_map;
            ///< An `LvaluePropertyMap` that maps to Face.
                
            using Element_with_properties = typename Element_map::key_type;
            ///< Value type of the iterator in the input range.

            /*!
                The constructor requires a mesh to perform neighbor search.
            */
            Connectivity_on_surface_mesh(const Mesh &mesh) :
            m_mesh(mesh) 
            { }

            /*!
                From a query element `face_index`, this function retrieves the face via the element map and uses `CGAL::faces_around_face` to get a list of neighbor faces. The result is returned in `neighbors`.
                \tparam Neighbors CGAL::Shape_detection::Region_growing::Neighbors
            */
            template<class Neighbors>
            void get_neighbors(const Element_with_properties &face_index, Neighbors &neighbors) {
                    
                Face face = get(m_elem_map, face_index);
                neighbors.clear();
                    
                const Face_range &local_faces = faces_around_face(m_mesh.halfedge(face), m_mesh);
                
                for (auto it = local_faces.begin(); it != local_faces.end(); ++it)
                    if ((*it) != Face(4294967295)) neighbors.push_back(*it);
            }

        private:

            // Fields.
            const Mesh &m_mesh;
            const Element_map m_elem_map = Element_map();
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_CONNECTIVITY_ON_SURFACE_MESH_H
