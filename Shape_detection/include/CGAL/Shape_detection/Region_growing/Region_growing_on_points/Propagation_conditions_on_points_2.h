#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_2_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_2_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include <CGAL/Shape_detection/Region_growing/Tools/Sqrt.h>
#include <CGAL/Shape_detection/Region_growing/Tools/Default_item_index_to_item_property_map.h>

namespace CGAL {

    namespace Shape_detection {

        /*! 
            \ingroup PkgShapeDetection
            \brief Local and global conditions for the region growing algorithm on a 2D point cloud.
            \tparam PointMap An `LvaluePropertyMap` that maps to `Point_2`.
            \tparam NormalMap An `LvaluePropertyMap` that maps to `Vector_2`.
            \tparam Traits Model of `RegionGrowingOnPointsTraits`
            \cgalModels `RegionGrowingPropagationConditions`
        */
        template<class InputRange, class PointMap, class NormalMap, class Traits, 
        class ItemIndexToItemMap = CGAL::Shape_detection::Default_item_index_to_item_property_map<InputRange> >
        class Propagation_conditions_on_points_2 {

        public:

            using Input_range            = InputRange;
            ///< An arbitrary range with user-defined items.

            using Point_map              = PointMap;
            ///< An `LvaluePropertyMap` that maps to `Point_2`.

            using Normal_map             = NormalMap;
            ///< An `LvaluePropertyMap` that maps to `Vector_2`.

            using Item_index_to_item_map = ItemIndexToItemMap;
            ///< An `LvaluePropertyMap` that maps to an arbitrary item.

            using FT                     = typename Traits::FT;       ///< Number type
            using Point_2                = typename Traits::Point_2;  ///< Point type
            using Vector_2               = typename Traits::Vector_2; ///< Vector type
            using Line_2                 = typename Traits::Line_2;   ///< Line type

            ///< \cond SKIP_IN_MANUAL
            using Local_traits           = Exact_predicates_inexact_constructions_kernel;
            using Local_FT               = typename Local_traits::FT;
            using Local_point_2          = typename Local_traits::Point_2;
            using Local_line_2           = typename Local_traits::Line_2;
            using To_local_converter     = Cartesian_converter<Traits, Local_traits>;

            using Squared_length_2       = typename Traits::Compute_squared_length_2;
            using Squared_distance_2     = typename Traits::Compute_squared_distance_2;
            using Scalar_product_2       = typename Traits::Compute_scalar_product_2;

            using Get_sqrt               = CGAL::Shape_detection::Get_sqrt<Traits>;
            using Sqrt                   = typename Get_sqrt::Sqrt;

            using Item_index             = std::size_t;
            ///< \endcond

            /*!
                Each region is represented by a line. The constructor requires three parameters, in order: 
                the maximum distance from a point to the region, 
                the minimum dot product between the normal associated with the point and the normal of the region, 
                and the minimum number of points a region must have.
                In addition, you can provide instances of the Point_map, Normal_map, and Traits classes.
            */
            Propagation_conditions_on_points_2(const Input_range &input_range, 
            const FT distance_threshold = FT(1), const FT normal_threshold = FT(9) / FT(10), const size_t min_region_size = 2, 
            const Point_map &point_map = PointMap(), const Normal_map &normal_map = Normal_map(), const Traits &traits = Traits()) : 
            m_item_index_to_item_map(input_range),
            m_distance_threshold(distance_threshold),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_point_map(point_map),
            m_normal_map(normal_map),
            m_squared_length_2(traits.compute_squared_length_2_object()),
            m_squared_distance_2(traits.compute_squared_distance_2_object()),
            m_scalar_product_2(traits.compute_scalar_product_2_object()),
            m_sqrt(Get_sqrt::sqrt_object(traits)) {

                CGAL_precondition(distance_threshold >= FT(0));
                CGAL_precondition(normal_threshold   >= FT(0) && normal_threshold <= FT(1));
                CGAL_precondition(min_region_size    > 1);
            }        

            Propagation_conditions_on_points_2(const Input_range &, const Item_index_to_item_map &item_index_to_item_map, 
            const FT distance_threshold = FT(1), const FT normal_threshold = FT(9) / FT(10), const size_t min_region_size = 2, 
            const Point_map &point_map = PointMap(), const Normal_map &normal_map = Normal_map(), const Traits &traits = Traits()) :
            m_item_index_to_item_map(item_index_to_item_map),
            m_distance_threshold(distance_threshold),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_point_map(point_map),
            m_normal_map(normal_map),
            m_squared_length_2(traits.compute_squared_length_2_object()),
            m_squared_distance_2(traits.compute_squared_distance_2_object()),
            m_scalar_product_2(traits.compute_scalar_product_2_object()),
            m_sqrt(Get_sqrt::sqrt_object(traits)) { 

                CGAL_precondition(distance_threshold >= FT(0));
                CGAL_precondition(normal_threshold   >= FT(0) && normal_threshold <= FT(1));
                CGAL_precondition(min_region_size    > 1);
            }

            /*!
                Local conditions that check if a query item belongs to the given region.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            bool belongs_to_region(const Item_index query_index, const Region &region) const {

                const auto &query_item = get(m_item_index_to_item_map, query_index);

                const Point_2  &query_point = get(m_point_map , *query_item);
                const Vector_2 &normal      = get(m_normal_map, *query_item);

                const FT normal_length = m_sqrt(m_squared_length_2(normal));
                CGAL_precondition(normal_length > FT(0));

                const Vector_2 query_normal = normal / normal_length;

                const FT distance_to_fitted_line = m_sqrt(m_squared_distance_2( query_point , m_line_of_best_fit));
                const FT cos_angle               = CGAL::abs(m_scalar_product_2(query_normal, m_normal_of_best_fit));

                return ( ( distance_to_fitted_line <= m_distance_threshold ) && ( cos_angle >= m_normal_threshold ) );
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
                Update the class's best fit line that will be used later by local conditions.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            void update(const Region &region) {

                CGAL_precondition(region.size() > 0);
                if (region.size() == 1) { // create new reference line and normal
                    
                    // The best fit line will be a line through this point with its normal being the point's normal.
                    const auto &item = get(m_item_index_to_item_map, region[0]);

                    const Point_2  &point  = get(m_point_map , *item);
                    const Vector_2 &normal = get(m_normal_map, *item);
                    
                    const FT normal_length = m_sqrt(m_squared_length_2(normal));
                    CGAL_precondition(normal_length > FT(0));

                    m_line_of_best_fit   = Line_2(point, normal).perpendicular(point);
                    m_normal_of_best_fit = normal / normal_length;

                } else { // update reference line and normal

                    Local_FT x = Local_FT(0);
                    Local_FT y = Local_FT(0);

                    std::vector<Local_point_2> points(region.size());
                    for (Item_index i = 0; i < region.size(); ++i) {

                        const auto &item = get(m_item_index_to_item_map, region[i]);
                        points[i] = m_to_local_converter(get(m_point_map, *item));

                        x += points[i].x();
                        y += points[i].y();
                    }

                    CGAL_precondition(points.size()> 0);

                    x /= static_cast<Local_FT>(points.size());
                    y /= static_cast<Local_FT>(points.size());

                    Local_line_2  fitted_line;
                    Local_point_2 fitted_centroid = Local_point_2(x, y);

                    // The best fit line will be a line fitted to all region points with its normal being perpendicular to the line.
                    #ifndef CGAL_EIGEN2_ENABLED
                        linear_least_squares_fitting_2(points.begin(), points.end(), fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 2>());
                    #else 
                        linear_least_squares_fitting_2(points.begin(), points.end(), fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 2>());
                    #endif
                    
                    m_line_of_best_fit = Line_2(static_cast<FT>(fitted_line.a()), static_cast<FT>(fitted_line.b()), static_cast<FT>(fitted_line.c()));
                    
                    const Vector_2 normal  = m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
                    const FT normal_length = m_sqrt(m_squared_length_2(normal));

                    CGAL_precondition(normal_length > FT(0));
                    m_normal_of_best_fit = normal / normal_length;
                }
            }

        private:
        
            // Fields.
            const FT                            m_distance_threshold;
            const FT                            m_normal_threshold;
            const size_t                        m_min_region_size;

            const Point_map                    &m_point_map;
            const Normal_map                   &m_normal_map;
            
            const Item_index_to_item_map        m_item_index_to_item_map;
            
            const Squared_length_2              m_squared_length_2;
            const Squared_distance_2            m_squared_distance_2;
            const Scalar_product_2              m_scalar_product_2;
            const Sqrt                          m_sqrt;

            const To_local_converter            m_to_local_converter;

            Line_2                              m_line_of_best_fit;
            Vector_2                            m_normal_of_best_fit;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_2_H
