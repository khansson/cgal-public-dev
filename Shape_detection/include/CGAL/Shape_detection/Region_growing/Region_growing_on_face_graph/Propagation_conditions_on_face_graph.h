#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_FACE_GRAPH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_FACE_GRAPH_H

// STL includes.
#include <vector>

// Boost includes.
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Face graph includes.
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Internal/Sqrt.h>
#include <CGAL/Shape_detection/Region_growing/Property_maps/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection.
            \brief Local and global conditions for the region growing algorithm on a face graph.
            \tparam FaceGraph General face graph. Model of `FaceGraph`.
            \tparam Mesh Has model `CGAL::Surface_mesh`.
            \cgalModels `RegionGrowingPropagationConditions`
        */
        template<class Traits, class FaceGraph,
        class FaceRange = typename FaceGraph::Face_range,
        class IndexToFaceMap   = CGAL::Shape_detection::Random_access_index_to_item_property_map<FaceRange>,
        class VertexToPointMap = typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type>
        class Propagation_conditions_on_face_graph {

        public:

            using Face_graph             = FaceGraph;
            ///< An arbitrary range with user-defined items.

            using Face_range             = FaceRange;

            using Index_to_face_map      = IndexToFaceMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

            using Vertex_to_point_map    = VertexToPointMap;

            using Face                   = typename Index_to_face_map::value_type;

            using FT                     = typename Traits::FT;       ///< Number type
            using Point_3                = typename Traits::Point_3;  ///< Point type
            using Vector_3               = typename Traits::Vector_3; ///< Vector type
            using Plane_3                = typename Traits::Plane_3;  ///< Plane type

            ///< \cond SKIP_IN_MANUAL
            using Local_traits           = Exact_predicates_inexact_constructions_kernel;
            using Local_FT               = typename Local_traits::FT;
            using Local_point_3          = typename Local_traits::Point_3;
            using Local_plane_3          = typename Local_traits::Plane_3;
            using To_local_converter     = Cartesian_converter<Traits, Local_traits>;

            using Squared_length_3       = typename Traits::Compute_squared_length_3;
            using Squared_distance_3     = typename Traits::Compute_squared_distance_3;
            using Scalar_product_3       = typename Traits::Compute_scalar_product_3;

            using Get_sqrt               = CGAL::Shape_detection::Get_sqrt<Traits>;
            using Sqrt                   = typename Get_sqrt::Sqrt;

            using Index                  = long;
            ///< \endcond

            /*!
                Each region is represented by a plane. The constructor requires three parameters, in order: 
                the maximum distance from a point to the region, 
                the minimum dot product between the normal associated with the point and the normal of the region, 
                and the minimum number of points a region must have.
                In addition, you can provide instances of the Face_descriptor_map and Traits classes.
            */
            Propagation_conditions_on_face_graph(const Face_graph &face_graph,
            const FT distance_threshold = FT(1), const FT normal_threshold = FT(9) / FT(10), const std::size_t min_region_size = 1, 
            const Vertex_to_point_map vertex_to_point_map = Vertex_to_point_map(), const Traits traits = Traits()) :
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(m_face_range),
            m_vertex_to_point_map(vertex_to_point_map),
            m_distance_threshold(distance_threshold),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_squared_length_3(traits.compute_squared_length_3_object()),
            m_squared_distance_3(traits.compute_squared_distance_3_object()),
            m_scalar_product_3(traits.compute_scalar_product_3_object()),
            m_sqrt(Get_sqrt::sqrt_object(traits)) {

                CGAL_precondition(distance_threshold >= FT(0));
                CGAL_precondition(normal_threshold   >= FT(0) && normal_threshold <= FT(1));
                CGAL_precondition(min_region_size    > 0);
            }        
            
            Propagation_conditions_on_face_graph(const Face_graph &face_graph, const Index_to_face_map index_to_face_map,
            const FT distance_threshold = FT(1), const FT normal_threshold = FT(9) / FT(10), const std::size_t min_region_size = 1, 
            const Vertex_to_point_map vertex_to_point_map = Vertex_to_point_map(), const Traits traits = Traits()) : 
            m_face_graph(face_graph),
            m_face_range(CGAL::faces(m_face_graph)),
            m_index_to_face_map(index_to_face_map),
            m_vertex_to_point_map(vertex_to_point_map),
            m_distance_threshold(distance_threshold),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_squared_length_3(traits.compute_squared_length_3_object()),
            m_squared_distance_3(traits.compute_squared_distance_3_object()),
            m_scalar_product_3(traits.compute_scalar_product_3_object()),
            m_sqrt(Get_sqrt::sqrt_object(traits)) {

                CGAL_precondition(distance_threshold >= FT(0));
                CGAL_precondition(normal_threshold   >= FT(0) && normal_threshold <= FT(1));
                CGAL_precondition(min_region_size    > 0);
            }    

            /*!
                Local conditions that check if a query item belongs to the given region.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            bool belongs_to_region(const Index query_index, const Region &region) const {

                const auto &face = get(m_index_to_face_map, query_index);

                Point_3 face_centroid;
                get_face_centroid(face, face_centroid);
                
                Vector_3 face_normal;
                get_face_normal(face, face_normal);

                const FT distance_to_fitted_plane = m_sqrt(m_squared_distance_3(face_centroid, m_plane_of_best_fit));
                const FT cos_angle                = CGAL::abs(m_scalar_product_3( face_normal, m_normal_of_best_fit));

                return ( ( distance_to_fitted_plane <= m_distance_threshold ) && ( cos_angle >= m_normal_threshold ) );
            }

            /*!
                Global conditions that check if a region size is large enough to be accepted.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            inline bool are_valid(const Region &region) const {
                return ( region.size() >= m_min_region_size );
            }

            /*!
                Update the class's best fit plane that will be used later by local conditions.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            void update(const Region &region) {

                CGAL_precondition(region.size() > 0);
                if (region.size() == 1) { // create new reference plane and normal
                        
                    const auto &face = get(m_index_to_face_map, region[0]);                    
                    Point_3 face_centroid;

                    get_face_centroid(face, face_centroid);
                    get_face_normal(face, m_normal_of_best_fit);

                    m_plane_of_best_fit  = Plane_3(face_centroid, m_normal_of_best_fit);

                } else {

                    Local_FT x = Local_FT(0);
                    Local_FT y = Local_FT(0);
                    Local_FT z = Local_FT(0);

                    std::vector<Local_point_3> points(region.size());
                    for (Index i = 0; i < region.size(); ++i) {
                        
                        const auto &face = get(m_index_to_face_map, region[i]);
                        const auto &halfedge = CGAL::halfedge(*face, m_face_graph);

                        const auto &vertices = CGAL::vertices_around_face(halfedge, m_face_graph);
                        for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex) {
                            
                            const Point_3 &tmp_point = get(m_vertex_to_point_map, *vertex);
                            points[i] = m_to_local_converter(tmp_point);

                            x += points[i].x();
                            y += points[i].y();
                            z += points[i].z();
                        }
                    }

                    CGAL_precondition(points.size() > 0);

                    x /= static_cast<Local_FT>(points.size());
                    y /= static_cast<Local_FT>(points.size());
                    z /= static_cast<Local_FT>(points.size());

                    Local_plane_3 fitted_plane;
                    Local_point_3 fitted_centroid = Local_point_3(x, y, z);

                    #ifndef CGAL_EIGEN3_ENABLED
                        linear_least_squares_fitting_3(points.begin(), points.end(), fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 3>());
                    #else 
                        linear_least_squares_fitting_3(points.begin(), points.end(), fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 3>());
                    #endif

                    m_plane_of_best_fit = Plane_3(static_cast<FT>(fitted_plane.a()), static_cast<FT>(fitted_plane.b()), static_cast<FT>(fitted_plane.c()), static_cast<FT>(fitted_plane.d()));

                    const Vector_3 normal  = m_plane_of_best_fit.orthogonal_vector();
                    const FT normal_length = m_sqrt(m_squared_length_3(normal));

                    CGAL_precondition(normal_length > FT(0));
                    m_normal_of_best_fit = normal / normal_length;
                }
            }

        private:
        
            void get_face_centroid(const Face &face, Point_3 &face_centroid) const {

                const auto &halfedge = CGAL::halfedge(*face, m_face_graph);
                const auto &vertices = CGAL::vertices_around_face(halfedge, m_face_graph);

                // Compute centroid.
                FT sum = FT(0);

                FT x = FT(0);
                FT y = FT(0);
                FT z = FT(0);
                
                for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex, sum += FT(1)) {
                    const Point_3 &point = get(m_vertex_to_point_map, *vertex);

                    x += point.x();
                    y += point.y();
                    z += point.z();
                }
                CGAL_precondition(sum > FT(0));

                x /= sum;
                y /= sum;
                z /= sum;

                face_centroid = Point_3(x, y, z);
            }

            void get_face_normal(const Face &face, Vector_3 &face_normal) const {

                const auto &halfedge = CGAL::halfedge(*face, m_face_graph);
                const auto &vertices = CGAL::vertices_around_face(halfedge, m_face_graph);

                auto vertex = vertices.begin();
                const Point_3 &p1 = get(m_vertex_to_point_map, *vertex); ++vertex;
                const Point_3 &p2 = get(m_vertex_to_point_map, *vertex); ++vertex;
                const Point_3 &p3 = get(m_vertex_to_point_map, *vertex);

                const Vector_3 tmp_normal  = CGAL::normal(p1, p2, p3);
                const FT tmp_normal_length = m_sqrt(m_squared_length_3(tmp_normal));

                CGAL_precondition(tmp_normal_length > FT(0));
                face_normal = tmp_normal / tmp_normal_length;
            }

            // Fields.
            const Face_graph                   &m_face_graph;
            const Face_range                    m_face_range;

            const Index_to_face_map             m_index_to_face_map;
            const Vertex_to_point_map           m_vertex_to_point_map;

            const FT                            m_distance_threshold;
            const FT                            m_normal_threshold;
            const std::size_t                   m_min_region_size;
            
            const Squared_length_3              m_squared_length_3;
            const Squared_distance_3            m_squared_distance_3;
            const Scalar_product_3              m_scalar_product_3;
            const Sqrt                          m_sqrt;

            const To_local_converter            m_to_local_converter;

            Plane_3                             m_plane_of_best_fit;
            Vector_3                            m_normal_of_best_fit;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_FACE_GRAPH_H
