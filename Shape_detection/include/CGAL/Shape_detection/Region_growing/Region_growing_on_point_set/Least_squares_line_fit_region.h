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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utilities.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*! 
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Best least squares line fit conditions on 2D points.

    This class implements propagation conditions for detecting lines 
    on 2D points via the `Shape_detection::Region_growing` approach, 
    where quality of detected lines is based on a local least squares line fit.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam InputRange 
    is a model of `ConstRange`. Its iterator type is `RandomAccessIterator`. 
    Its value type depends on the item type used in Region Growing, 
    for example it can be `std::pair<CGAL::Point_2, int>` 
    or any user-defined type.

    \tparam PointMap 
    is an `LvaluePropertyMap` that maps to `CGAL::Point_2`.

    \tparam NormalMap 
    is an `LvaluePropertyMap` that maps to `CGAL::Vector_2`.
    
    \cgalModels `RegionGrowingPropagationConditions`
  */
  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap, 
  typename NormalMap>
  class Least_squares_line_fit_region {

  public:

    /// \name Types
    /// @{
    
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Normal_map = NormalMap;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// Vector type.
    typedef typename GeomTraits::Vector_2 Vector_2;

    /// Type of the line.
    typedef typename GeomTraits::Line_2 Line_2;

    /// \cond SKIP_IN_MANUAL
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_2 = typename Local_traits::Point_2;
    using Local_line_2 = typename Local_traits::Line_2;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_2 = typename Traits::Compute_squared_length_2;
    using Squared_distance_2 = typename Traits::Compute_squared_distance_2;
    using Scalar_product_2 = typename Traits::Compute_scalar_product_2;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param input_range 
      An instance of an `InputRange` container with 2D points and 
      2D associated normal vectors.

      \param distance_threshold 
      Maximum distance from a point to the region represented by a line of 
      the type `CGAL::Line_2`.

      \param normal_threshold 
      Minimum dot product between the normal assigned to the point and 
      the normal assigned to the region represented by a line 
      of the type `CGAL::Line_2`.

      \param min_region_size 
      The minimum number of 2D points a region must have.

      \param point_map
      An instance of an `LvaluePropertyMap` that maps an item from `input_range` 
      to `CGAL::Point_2`.

      \param normal_map
      An instance of an `LvaluePropertyMap` that maps an item from `input_range` 
      to `CGAL::Vector_2`.

      \param traits
      An instance of the `GeomTraits` class.

      \pre `input_range.size() > 0`
      \pre `distance_threshold >= 0`
      \pre `normal_threshold >= 0 && normal_threshold <= 1`
      \pre `min_region_size > 0`
    */
    Least_squares_line_fit_region(
      const InputRange& input_range, 
      const FT distance_threshold = FT(1), 
      const FT angle_threshold = FT(25), 
      const std::size_t min_region_size = 2, 
      const PointMap point_map = PointMap(), 
      const NormalMap normal_map = NormalMap(), 
      const GeomTraits traits = GeomTraits()) :
    m_input_range(input_range),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_min_region_size(min_region_size),
    m_point_map(point_map),
    m_normal_map(normal_map),
    m_squared_length_2(traits.compute_squared_length_2_object()),
    m_squared_distance_2(traits.compute_squared_distance_2_object()),
    m_scalar_product_2(traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(input_range.size() > 0); 

      CGAL_precondition(distance_threshold >= FT(0));
      CGAL_precondition(angle_threshold >= FT(0) && angle_threshold <= FT(90));
      CGAL_precondition(min_region_size > 0);
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief Checks if a point belongs to a region.

      Checks if the point with the index `query_index` belongs to the region
      that is currently getting developed using the `distance_threshold` and
      `normal_threshold` values.

      \param query_index
      Index of the query point.

      The second parameter is not used in this implementation.

      Implements the function `RegionGrowingPropagationConditions::belongs_to_region()`.

      \pre `query_index >= 0 && query_index < total_number_of_points`
    */
    bool is_part_of_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      const auto& key = *(m_input_range.begin() + query_index);
      const Point_2& query_point = get(m_point_map, key);

      const Vector_2& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_2(normal));
      CGAL_precondition(normal_length > FT(0));
      const Vector_2 query_normal = normal / normal_length;

      const FT distance_to_fitted_line = 
      m_sqrt(m_squared_distance_2(query_point, m_line_of_best_fit));
      
      const FT cos_value = 
      CGAL::abs(m_scalar_product_2(query_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_line <= m_distance_threshold ) && 
        ( cos_value >= m_normal_threshold ));
    }

    /*!
      \brief Validates the final `region`.

      Controls if the `region` that has been created contains at least 
      `min_region_size` points.

      \param region
      Stores indices of all points that belong to the region.

      Implements the function `RegionGrowingPropagationConditions::is_valid_region()`.
    */
    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return ( region.size() >= m_min_region_size );
    }

    /*!
      \brief Recomputes the least squares line.

      Recomputes the internal least squares line that represents the `region`
      currently being developed. The line is fitted to all points, 
      which have been added to the `region` so far.

      \param region
      Stores indices of all points that belong to the region.

      Implements the function `RegionGrowingPropagationConditions::update()`.

      \pre `region.size() > 0`
    */
    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and normal
                    
        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());

        // The best fit line will be a line through this point with 
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + region[0]);

        const Point_2& point = get(m_point_map, key);
        const Vector_2& normal = get(m_normal_map, key);
                    
        const FT normal_length = m_sqrt(m_squared_length_2(normal));
        CGAL_precondition(normal_length > FT(0));

        m_normal_of_best_fit = 
        normal / normal_length;
        
        m_line_of_best_fit = 
        Line_2(point, m_normal_of_best_fit).perpendicular(point);

      } else { // update reference line and normal

        std::vector<Local_point_2> points;
        points.reserve(region.size());

        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + region[i]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == region.size());

        Local_line_2 fitted_line;
        Local_point_2 fitted_centroid;

        // The best fit line will be a line fitted to all region points with 
        // its normal being perpendicular to the line.
        #ifndef CGAL_EIGEN2_ENABLED
          linear_least_squares_fitting_2(
            points.begin(), points.end(), 
            fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 2>());
        #else 
          linear_least_squares_fitting_2(
            points.begin(), points.end(), 
            fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 2>());
        #endif
                    
        m_line_of_best_fit = 
        Line_2(
          static_cast<FT>(fitted_line.a()), 
          static_cast<FT>(fitted_line.b()), 
          static_cast<FT>(fitted_line.c()));
                    
        const Vector_2 normal = 
        m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
        const FT normal_length = m_sqrt(m_squared_length_2(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

    /// @}

  private:
        
    // Fields.
    const Input_range& m_input_range;
            
    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const std::size_t m_min_region_size;

    const Point_map m_point_map;
    const Normal_map m_normal_map;
            
    const Squared_length_2 m_squared_length_2;
    const Squared_distance_2 m_squared_distance_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Line_2 m_line_of_best_fit;
    Vector_2 m_normal_of_best_fit;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
