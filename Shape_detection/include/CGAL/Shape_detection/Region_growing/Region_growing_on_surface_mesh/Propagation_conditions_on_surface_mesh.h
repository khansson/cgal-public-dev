#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_SURFACE_MESH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_SURFACE_MESH_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Todo:
// Add default values to the constructor. Also specify preconditions on these parameters.
// Change global CGAL functions like squared_distance to their counterparts from the Kernel.
// Change m_line_of_best_fit and m_normal_of_best_fit to the FT type? - Does it actually work with exact Kernel?
// Add m_min_region_size parameter.
// Add const to Vertex_range vr.
// Add const to normal.
// Add const and reference to Point_3 tmp.
// Remove conversion between vector and point in the find_centroid() function.

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection.
            \brief Local and global conditions for the region growing algorithm on a mesh.
            \tparam RegionGrowingTraits `CGAL::Shape_detection::Region_growing_traits`
            \tparam Mesh Has model `CGAL::Surface_mesh`.
            \cgalModels `RegionGrowingPropagationConditions`
        */
        template<class RegionGrowingTraits, class Mesh>
        class Propagation_conditions_on_surface_mesh {

            // Sqrt object.
            template<class ...>
            using void_t = void;

            template<class Kernel, class = void>
            class Get_sqrt {
                using FT = typename Kernel::FT;
            public:
                FT operator()(const FT value) const {
                    return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
                }
            };

            template<class Kernel>
            class Get_sqrt<Kernel, void_t<typename Kernel::Sqrt> > : 
            Kernel::Sqrt 
            { };

        public:

            using Kernel = typename RegionGrowingTraits::Kernel;
            
            using Vertex_index = typename Mesh::Vertex_index;

            using Vertex_range = Iterator_range< Vertex_around_face_iterator<Mesh> >;
            ///< Result type of CGAL::vertices_around_face().
            
            using Face = typename RegionGrowingTraits::Element;
            ///< Must be equivalent to Mesh::Face_index.

            using Element_map = typename RegionGrowingTraits::Element_map; 
            ///< An `LvaluePropertyMap` that maps to Face.

            using Element_with_properties = typename Element_map::key_type;
            ///< Value type of the iterator in the input range.

            using FT       = typename Kernel::FT;       ///< Number type
            using Point_3  = typename Kernel::Point_3;  ///< Point type
            using Vector_3 = typename Kernel::Vector_3; ///< Vector type
            using Plane_3  = typename Kernel::Plane_3;  ///< Plane type

            ///< \cond SKIP_IN_MANUAL
            using Sqrt               = Get_sqrt<Kernel>;
            using Local_kernel       = Exact_predicates_inexact_constructions_kernel;
            using Local_FT           = Local_kernel::FT;
            using Local_point_3      = Local_kernel::Point_3;
            using Local_vector_3     = Local_kernel::Vector_3;
            using Local_plane_3      = Local_kernel::Plane_3;
            using To_local_converter = Cartesian_converter<Kernel, Local_kernel>;
            ///< \endcond

            /*!
                Each region is represented by a plane. The constructor requires three parameters, in order: the mesh, the maximum distance from a point to the region, and the minimum dot product between the normal associated with the point and the normal of the region.
            */                
            Propagation_conditions_on_surface_mesh(const Mesh &mesh, const FT epsilon, const FT normal_threshold) :
            m_mesh(mesh),
            m_epsilon(epsilon),
            m_normal_threshold(normal_threshold) 
            { }
                
            Propagation_conditions_on_surface_mesh(const Mesh &mesh, const FT epsilon, const FT normal_threshold,
            const Element_map &element_map) :
            m_mesh(mesh),
            m_epsilon(epsilon),
            m_normal_threshold(normal_threshold),
            m_element_map(element_map)
            { }

            /*!
                Local conditions that check if a new face in `unassigned_element` is similar to the face `assigned_element` and its enclosing region `region`.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            bool are_in_same_region(
                                const Element_with_properties &assigned_element,
                                const Element_with_properties &unassigned_element,
                                const Region &region) {

                const Face &face = get(m_element_map, unassigned_element);

                Point_3 face_centroid;
                get_centroid(face, face_centroid);
                
                Vector_3 face_normal;
                get_normal(face, face_normal);

                // Must use Local_FT, because fit plane is of local kernel.
                const Local_FT distance_to_fit_plane = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(face_centroid), m_plane_of_best_fit));
                const Local_FT cos_angle             = CGAL::abs(m_to_local_converter(face_normal) * m_normal_of_best_fit);

                CGAL_precondition(m_epsilon >= FT(0));
                CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));

                return ( ( distance_to_fit_plane <= m_to_local_converter(m_epsilon) ) && ( cos_angle >= m_to_local_converter(m_normal_threshold) ) );
            }

            /*!
                Global conditions that check the validity of a region. Always return true if the region is not empty.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            inline bool are_valid(const Region &region) {
                
                return region.size() > 0;
            }

            /*!
                Update the class's best fit plane that will be used later by local conditions.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            void update(const Region &region) {

                CGAL_precondition(region.end() != region.begin());
                if (region.size() == 1) {
                        
                    const Face &face = get(m_element_map, *region.begin());

                    Point_3 face_centroid;
                    get_centroid(face, face_centroid);
                    
                    Vector_3 face_normal;
                    get_normal(face, face_normal);

                    const FT normal_length = m_sqrt(face_normal.squared_length());

                    m_plane_of_best_fit  = m_to_local_converter(Plane_3(face_centroid, face_normal));
                    m_normal_of_best_fit = m_to_local_converter(face_normal / normal_length);

                } else {

                    // Extract the points of the region (Point_3).
                    std::vector<Local_point_3> points;

                    // The region is formed by face indices.
                    for (auto rit = region.begin(); rit != region.end(); ++rit) {
                        const Face &face = get(m_element_map, *rit);
                            
                        // Get vertices of each face and push them to `points`.
                        Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);

                        for (auto vit = vr.begin(); vit != vr.end(); ++vit)
                            points.push_back(m_to_local_converter(m_mesh.point(*vit)));
                    }
                        
                    // Fit the region to a plane.
                    Local_point_3 centroid;

                    #ifndef CGAL_EIGEN3_ENABLED
                        linear_least_squares_fitting_3(points.begin(), points.end(), m_plane_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), CGAL::Default_diagonalize_traits<Local_FT, 3>());
                    #else 
                        linear_least_squares_fitting_3(points.begin(), points.end(), m_plane_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), CGAL::Eigen_diagonalize_traits<Local_FT, 3>());
                    #endif

                    Local_vector_3 normal = m_plane_of_best_fit.orthogonal_vector();

                    const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                    m_normal_of_best_fit = normal / normal_length;
                }
            }

        private:

            // Functions.
            void get_centroid(const Face &face, Point_3 &centroid) const {
                Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);

                // Calculate centroid.
                size_t n = 0;
                Vector_3 center = Vector_3(0, 0, 0); // temporarily set as a vector to perform + and / operators
                
                for (auto it = vr.begin(); it != vr.end(); ++it, ++n) {

                    Point_3 point = m_mesh.point(*it);
                    center += Vector_3(point.x(), point.y(), point.z());
                }
                
                CGAL_precondition(n > 2);
                center /= n;

                centroid = Point_3(center.x(), center.y(), center.z()); // convert back to point
            }

            void get_normal(const Face &face, Vector_3 &normal) const {

                Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);

                auto it = vr.begin();
                normal = CGAL::normal(m_mesh.point(*it++), m_mesh.point(*it++), m_mesh.point(*it++));
            }

            // Fields.
            const Mesh                     &m_mesh;
            
            const FT                        m_epsilon;
            const FT                        m_normal_threshold;

            const Sqrt                      m_sqrt;
            const Element_map               m_element_map;
            const To_local_converter        m_to_local_converter;
            
            Local_plane_3                   m_plane_of_best_fit;
            Local_vector_3                  m_normal_of_best_fit;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_SURFACE_MESH_H
