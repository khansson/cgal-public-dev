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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utilities.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*! 
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Least squares plane fit conditions on 3D points.

    This class implements propagation conditions for detecting planes 
    on 3D points via the `Shape_detection::Region_growing` approach, 
    where quality of detected planes is based on a local least squares plane fit.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam InputRange 
    is a model of `ConstRange`. Its iterator type is `RandomAccessIterator`. 
    Its value type depends on the item type used in Region Growing, 
    for example it can be `std::pair<CGAL::Point_3, CGAL::Vector_3>` 
    or any user-defined type.

    \tparam PointMap 
    is an `LvaluePropertyMap` that maps key to `CGAL::Point_3`.

    \tparam NormalMap 
    is an `LvaluePropertyMap` that maps to `CGAL::Vector_3`.
    
    \cgalModels `RegionGrowingPropagationConditions`
  */
  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap, 
  typename NormalMap>
  class Least_squares_plane_fit_region {

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
    typedef typename GeomTraits::Point_3 Point_3;

    /// Vector type.
    typedef typename GeomTraits::Vector_3 Vector_3;

    /// Type of the plane.
    typedef typename GeomTraits::Plane_3 Plane_3;

    /// \cond SKIP_IN_MANUAL
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_3 = typename Traits::Compute_squared_length_3;
    using Squared_distance_3 = typename Traits::Compute_squared_distance_3;
    using Scalar_product_3 = typename Traits::Compute_scalar_product_3;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param input_range 
      An instance of an `InputRange` container with 3D points 
      and 3D associated normal vectors.

      \param distance_threshold 
      Maximum distance from a point to the region represented by a plane of 
      the type `CGAL::Plane_3`.

      \param angle_threshold 
      Minimum dot product between the normal assigned to the point and 
      the normal assigned to the region represented by a plane 
      of the type `CGAL::Plane_3`.

      \param min_region_size 
      The minimum number of 3D points a region must have.

      \param point_map
      An instance of an `LvaluePropertyMap` that maps an item from `input_range` 
      to `CGAL::Point_3`.

      \param normal_map
      An instance of an `LvaluePropertyMap` that maps an item from `input_range` 
      to `CGAL::Vector_3`.

      \param traits
      An instance of the `GeomTraits` class.

      \pre `input_range.size() > 0`
      \pre `distance_threshold >= 0`
      \pre `normal_threshold >= 0 && normal_threshold <= 1`
      \pre `min_region_size > 0`
    */
    Least_squares_plane_fit_region(
      const InputRange& input_range, 
      const FT distance_threshold = FT(1), 
      const FT angle_threshold = FT(25), 
      const std::size_t min_region_size = 3, 
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
    m_squared_length_3(traits.compute_squared_length_3_object()),
    m_squared_distance_3(traits.compute_squared_distance_3_object()),
    m_scalar_product_3(traits.compute_scalar_product_3_object()),
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
      const Point_3& query_point = get(m_point_map, key);
      
      const Vector_3& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_3(normal));
      CGAL_precondition(normal_length > FT(0));
      const Vector_3 query_normal = normal / normal_length;

      const FT distance_to_fitted_plane = 
      m_sqrt(m_squared_distance_3(query_point, m_plane_of_best_fit));
      
      const FT cos_value = 
      CGAL::abs(m_scalar_product_3(query_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_plane <= m_distance_threshold ) && 
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
      \brief Recomputes the least squares plane.

      Recomputes the internal least squares plane that represents the `region`
      currently being developed. The plane is fitted to all points, 
      which have been added to the `region` so far.

      \param region
      Stores indices of all points that belong to the region.

      Implements the function `RegionGrowingPropagationConditions::update()`.

      \pre `region.size() > 0`
    */
    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference plane and normal
                    
        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());

        // The best fit plane will be a plane through this point with 
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + region[0]);

        const Point_3& point = get(m_point_map, key);
        const Vector_3& normal = get(m_normal_map, key);
                    
        const FT normal_length = m_sqrt(m_squared_length_3(normal));
        CGAL_precondition(normal_length > FT(0));

        m_normal_of_best_fit = 
        normal / normal_length;
        
        m_plane_of_best_fit = 
        Plane_3(point, m_normal_of_best_fit);

      } else { // update reference plane and normal

        std::vector<Local_point_3> points;
        points.reserve(region.size());

        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + region[i]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == region.size());

        Local_plane_3 fitted_plane;
        Local_point_3 fitted_centroid;

        // The best fit plane will be a plane fitted to all region points with 
        // its normal being perpendicular to the plane.
        CGAL::linear_least_squares_fitting_3(
          points.begin(), points.end(), 
          fitted_plane, fitted_centroid, 
          CGAL::Dimension_tag<0>());
                    
        m_plane_of_best_fit = 
        Plane_3(
          static_cast<FT>(fitted_plane.a()), 
          static_cast<FT>(fitted_plane.b()), 
          static_cast<FT>(fitted_plane.c()), 
          static_cast<FT>(fitted_plane.d()));
                    
        const Vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
        const FT normal_length = m_sqrt(m_squared_length_3(normal));

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
            
    const Squared_length_3 m_squared_length_3;
    const Squared_distance_3 m_squared_distance_3;
    const Scalar_product_3 m_scalar_product_3;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Plane_3 m_plane_of_best_fit;
    Vector_3 m_normal_of_best_fit;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_REGION_H
