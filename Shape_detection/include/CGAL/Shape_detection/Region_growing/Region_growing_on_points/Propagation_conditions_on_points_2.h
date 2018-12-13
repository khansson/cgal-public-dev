#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_2_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_2_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Todo:
// Add default values to the constructor. Also specify preconditions on these parameters.
// Change global CGAL functions like squared_distance to their counterparts from the Kernel.
// Change m_line_of_best_fit and m_normal_of_best_fit to the FT type? - Does it actually work with exact Kernel?
// Add const to normal_unassigned.
// Add const to normal.

namespace CGAL {

    namespace Shape_detection {

        /*! 
            \ingroup PkgShapeDetection
            \brief Local and global conditions for the region growing algorithm on a 2D point cloud.
            \tparam RegionGrowingTraits CGAL::Shape_detection::Region_growing_traits
            \tparam NormalMap An `LvaluePropertyMap` that maps to a vector, representing the normal associated with the point.
            \cgalModels `RegionGrowingPropagationConditions`
        */
        template<class RegionGrowingTraits, class NormalMap>
        class Propagation_conditions_on_points_2 {

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
            
            using Normal_map = NormalMap;
            
            using Element_map = typename RegionGrowingTraits::Element_map;

            using Element_with_properties = typename Element_map::key_type;
            ///< Value type of the iterator in the input range.

            using FT       = typename Kernel::FT;       ///< Number type
            using Point_2  = typename Kernel::Point_2;  ///< Point type
            using Vector_2 = typename Kernel::Vector_2; ///< Vector type
            using Line_2   = typename Kernel::Line_2;   ///< Line type

            ///< \cond SKIP_IN_MANUAL
            using Sqrt               = Get_sqrt<Kernel>;
            using Local_kernel       = Exact_predicates_inexact_constructions_kernel;
            using Local_FT           = typename Local_kernel::FT;
            using Local_point_2      = typename Local_kernel::Point_2;
            using Local_vector_2     = typename Local_kernel::Vector_2;
            using Local_line_2       = typename Local_kernel::Line_2;
            using To_local_converter = Cartesian_converter<Kernel, Local_kernel>;
            ///< \endcond

            /*!
                Each region is represented by a line. The constructor requires three parameters, in order: the maximum distance from a point to the region, the minimum dot product between the normal associated with the point and the normal of the region, and the minimum number of points a region must have.
            */
            Propagation_conditions_on_points_2(const FT epsilon, const FT normal_threshold, const size_t min_region_size) :
            m_epsilon(epsilon),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_sqrt(Sqrt()) 
            { }

            Propagation_conditions_on_points_2(const FT epsilon, const FT normal_threshold, const size_t min_region_size,
            const Element_map &element_map, const Normal_map &normal_map) :
            m_epsilon(epsilon),
            m_normal_threshold(normal_threshold),
            m_min_region_size(min_region_size),
            m_sqrt(Sqrt()),
            m_element_map(element_map),
            m_normal_map(normal_map)
            { }

            /*!
                Local conditions that check if a new point in `unassigned_element` is similar to the point `assigned_element` and its enclosing region. Each item in `Region` has the same structure as in the input range, i.e. has the type `Element_with_properties`.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Input_element, class Region>
            bool are_in_same_region(
                                const Input_element &assigned_element,
                                const Input_element &unassigned_element,
                                const Region &region) {

                const Point_2  &point_unassigned = get(m_element_map, *unassigned_element);
                const Vector_2 &normal           = get(m_normal_map , *unassigned_element);

                const FT normal_length     = m_sqrt(normal.squared_length());
                Vector_2 normal_unassigned = normal / normal_length;

                // Must use Local_FT, because fit line is of local kernel.
                const Local_FT distance_to_fit_line = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(point_unassigned), m_line_of_best_fit));
                const Local_FT cos_angle            = CGAL::abs(m_to_local_converter(normal_unassigned) * m_normal_of_best_fit);

                CGAL_precondition(m_epsilon >= FT(0));
                CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));

                return ( ( distance_to_fit_line <= m_to_local_converter(m_epsilon) ) && ( cos_angle >= m_to_local_converter(m_normal_threshold) ) );
            }

            /*!
                Global conditions that check if a region size is large enough to be accepted.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            inline bool are_valid(const Region &region) const {

                CGAL_precondition(m_min_region_size >= 2);
                return ( region.size() >= m_min_region_size );
            }

            /*!
                Update the class's best fit line that will be used later by local conditions.
                \tparam Region CGAL::Shape_detection::Region_growing::Region
            */
            template<class Region>
            void update(const Region &region) {

                CGAL_precondition(region.end() != region.begin());
                if (region.size() == 1) {
                    
                    // The best fit line will be a line through this point with its normal being the point's normal.
                    const Point_2  &point  = get(m_element_map, **region.begin());
                    const Vector_2 &normal = get(m_normal_map , **region.begin());
                    
                    const FT normal_length = m_sqrt(normal.squared_length());

                    m_line_of_best_fit   = m_to_local_converter(Line_2(point, normal).perpendicular(point));
                    m_normal_of_best_fit = m_to_local_converter(normal / normal_length);

                } else {

                    // Extract the geometric Element (Point_2).
                    size_t i = 0;
                    std::vector<Local_point_2> points(region.size());

                    for (auto it = region.begin(); it != region.end(); ++it, ++i)
                        points[i] = m_to_local_converter(get(m_element_map, **it));

                    // Fit the region to a line.
                    Local_point_2 centroid;

                    #ifndef CGAL_EIGEN2_ENABLED
                        linear_least_squares_fitting_2(points.begin(), points.end(), m_line_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), CGAL::Default_diagonalize_traits<Local_FT, 2>());
                    #else 
                        linear_least_squares_fitting_2(points.begin(), points.end(), m_line_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), CGAL::Eigen_diagonalize_traits<Local_FT, 2>());
                    #endif

                    Local_vector_2 normal = m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
                    
                    const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                    m_normal_of_best_fit = normal / normal_length;
                }
            }

        private:
        
            // Fields.
            const FT                        m_epsilon;
            const FT                        m_normal_threshold;
            const size_t                    m_min_region_size;

            const Sqrt                      m_sqrt;

            const Element_map               m_element_map;
            const Normal_map                m_normal_map;

            const To_local_converter        m_to_local_converter;

            Local_line_2                    m_line_of_best_fit;
            Local_vector_2                  m_normal_of_best_fit;
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_PROPAGATION_CONDITIONS_ON_POINTS_2_H
