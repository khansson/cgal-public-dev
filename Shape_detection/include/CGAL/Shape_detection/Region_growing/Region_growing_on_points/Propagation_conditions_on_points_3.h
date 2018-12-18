#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_3_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Tools/Sqrt.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Random_access_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*! 
            \ingroup PkgShapeDetection
            \brief Local and global conditions for the region growing algorithm on a 3D point cloud.
            \tparam PointMap An `LvaluePropertyMap` that maps to `Point_3`.
            \tparam NormalMap An `LvaluePropertyMap` that maps to `Vector_3`.
            \tparam Traits Model of `RegionGrowingOnPointsTraits`
            \cgalModels `RegionGrowingPropagationConditions`
        */
        template<class InputRange, class PointMap, class NormalMap, class Traits, 
        class IndexToItemMap = CGAL::Shape_detection::Random_access_index_to_item_property_map<InputRange> >
        class Propagation_conditions_on_points_3 {

        public:

            using Input_range            = InputRange;
            ///< An arbitrary range with user-defined items.

            using Point_map              = PointMap;
            ///< An `LvaluePropertyMap` that maps to `Point_3`.

            using Normal_map             = NormalMap;
            ///< An `LvaluePropertyMap` that maps to `Vector_3`.

            using Index_to_item_map      = IndexToItemMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

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

            using Index                  = std::size_t;
            ///< \endcond

            /*!
                Each region is represented by a plane. The constructor requires three parameters, in order: 
                the maximum distance from a point to the region, 
                the minimum dot product between the normal associated with the point and the normal of the region, 
                and the minimum number of points a region must have.
                In addition, you can provide instances of the Point_map, Normal_map, and Traits classes.
            */
            Propagation_conditions_on_points_3(const Input_range &input_range, 
            const FT distance_threshold = FT(1), const FT normal_threshold = FT(9) / FT(10), const std::size_t min_region_size = 3, 
            const Point_map point_map = Point_map(), const Normal_map normal_map = Normal_map(), const Traits traits = Traits()) : 
            m_input_range(input_range),
            m_index_to_item_map(m_input_range),
            m_distance_threshold(distance_threshold),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_point_map(point_map),
            m_normal_map(normal_map),
            m_squared_length_3(traits.compute_squared_length_3_object()),
            m_squared_distance_3(traits.compute_squared_distance_3_object()),
            m_scalar_product_3(traits.compute_scalar_product_3_object()),
            m_sqrt(Get_sqrt::sqrt_object(traits)) {

                CGAL_precondition(distance_threshold >= FT(0));
                CGAL_precondition(normal_threshold   >= FT(0) && normal_threshold <= FT(1));
                CGAL_precondition(min_region_size    > 2);
            }        

            Propagation_conditions_on_points_3(const Input_range &input_range, const Index_to_item_map index_to_item_map,
            const FT distance_threshold = FT(1), const FT normal_threshold = FT(9) / FT(10), const std::size_t min_region_size = 3, 
            const Point_map point_map = Point_map(), const Normal_map normal_map = Normal_map(), const Traits traits = Traits()) :
            m_input_range(input_range),
            m_index_to_item_map(index_to_item_map),
            m_distance_threshold(distance_threshold),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_point_map(point_map),
            m_normal_map(normal_map),
            m_squared_length_3(traits.compute_squared_length_3_object()),
            m_squared_distance_3(traits.compute_squared_distance_3_object()),
            m_scalar_product_3(traits.compute_scalar_product_3_object()),
            m_sqrt(Get_sqrt::sqrt_object(traits)) { 

                CGAL_precondition(distance_threshold >= FT(0));
                CGAL_precondition(normal_threshold   >= FT(0) && normal_threshold <= FT(1));
                CGAL_precondition(min_region_size    > 2);
            }

            /*!
                Local conditions that check if a query item belongs to the given region.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            bool belongs_to_region(const Index query_index, const Region &region) const {

                const auto &query_item = get(m_index_to_item_map, query_index);
                const auto &key        = *query_item;

                const Point_3  &query_point = get(m_point_map , key);
                const Vector_3 &normal      = get(m_normal_map, key);

                const FT normal_length = m_sqrt(m_squared_length_3(normal));
                CGAL_precondition(normal_length > FT(0));

                const Vector_3 query_normal = normal / normal_length;

                const FT distance_to_fitted_plane = m_sqrt(m_squared_distance_3( query_point , m_plane_of_best_fit));
                const FT cos_angle                = CGAL::abs(m_scalar_product_3(query_normal, m_normal_of_best_fit));

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
                    
                    // The best fit plane will be a plane through this point with its normal being the point's normal.
                    const auto &item = get(m_index_to_item_map, region[0]);
                    const auto &key  = *item;

                    const Point_3  &point  = get(m_point_map , key);
                    const Vector_3 &normal = get(m_normal_map, key);
                    
                    const FT normal_length = m_sqrt(m_squared_length_3(normal));
                    CGAL_precondition(normal_length > FT(0));

                    m_plane_of_best_fit  = Plane_3(point, normal);
                    m_normal_of_best_fit = normal / normal_length;

                } else { // update reference plane and normal

                    Local_FT x = Local_FT(0);
                    Local_FT y = Local_FT(0);
                    Local_FT z = Local_FT(0);

                    std::vector<Local_point_3> points(region.size());
                    for (Index i = 0; i < region.size(); ++i) {

                        const auto &item = get(m_index_to_item_map, region[i]);
                        const auto &key  = *item;

                        points[i] = m_to_local_converter(get(m_point_map, key));

                        x += points[i].x();
                        y += points[i].y();
                        z += points[i].z();
                    }

                    CGAL_precondition(points.size() > 0);

                    x /= static_cast<Local_FT>(points.size());
                    y /= static_cast<Local_FT>(points.size());
                    z /= static_cast<Local_FT>(points.size());

                    Local_plane_3 fitted_plane;
                    Local_point_3 fitted_centroid = Local_point_3(x, y, z);

                    // The best fit plane will be a plane fitted to all region points with its normal being perpendicular to the plane.
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
        
            // Fields.
            const Input_range                  &m_input_range;

            const Index_to_item_map             m_index_to_item_map;

            const FT                            m_distance_threshold;
            const FT                            m_normal_threshold;
            const std::size_t                   m_min_region_size;

            const Point_map                     m_point_map;
            const Normal_map                    m_normal_map;
            
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

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_3_H
