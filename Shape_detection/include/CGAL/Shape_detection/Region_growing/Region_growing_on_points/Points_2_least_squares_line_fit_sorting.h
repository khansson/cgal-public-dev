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

// #include <CGAL/license/Shape_detection.h>

// STL includes.
#include <map>
#include <vector>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_2_LEAST_SQUARES_LINE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINTS_2_LEAST_SQUARES_LINE_FIT_SORTING_H

namespace CGAL {
namespace Shape_detection {

  template<
    class GeomTraits,
    class InputRange,
    class PointMap,
    class Connectivity_>
  class Points_2_least_squares_line_fit_sorting {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Input_connectivity = Connectivity_;

    class Seed_map {
                        
    public:
      using key_type = std::size_t;
      using value_type = std::size_t;
      using category = boost::lvalue_property_map_tag;

      Seed_map(const std::vector<std::size_t>& objects_map) : m_objects_map(objects_map) { }

      value_type operator[](const key_type key) const { 
        return m_objects_map[key];
      }

      friend value_type get(
        const Seed_map& seed_map, 
        const key_type key) { 
        
        return seed_map[key];
      }

    private:
      const std::vector<std::size_t>& m_objects_map;
    };

    Points_2_least_squares_line_fit_sorting(
      const Input_range& input_range,
      Input_connectivity& connectivity,
      const Point_map point_map = Point_map()) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_connectivity(connectivity) { 
      std::size_t input_size = m_input_range.end() - m_input_range.begin();
      m_order.resize(input_size);
      m_scores.resize(input_size);
      for (std::size_t i = 0; i < input_size; ++i) m_order[i] = i;
    }

    void sort() {
      calculate_scores();
      Compare cmp(m_scores);
      // std::cerr << "* sorting\n";
      std::sort(m_order.begin(), m_order.end(), cmp);
      // std::cerr << "* sorted\n";
      // m_seed_map = Seed_map(m_order);
    }

    Seed_map seed_map() {
      return Seed_map(m_order);
    }

  private:

    /// \cond SKIP_IN_MANUAL
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_2 = typename Local_traits::Point_2;
    using Local_line_2 = typename Local_traits::Line_2;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;
    /// \endcond

    struct Compare {

      Compare(const std::vector<Local_FT>& scores) : m_scores(scores) { }

      bool operator()(const std::size_t x, const std::size_t y) const {
        // x stands before y in m_order iff m_scores[x] > m_scores[y]
        return m_scores[x] > m_scores[y];
      }

      private:
        const std::vector<Local_FT>& m_scores;

    };

    void calculate_scores() {
      // std::cerr << "* calculating scores\n";
      // std::cerr << m_scores.size() << std::endl;
      std::size_t input_size = m_input_range.end() - m_input_range.begin();
      for (int i = 0; i < input_size; ++i) {
        std::vector<std::size_t> neighbors;
        m_connectivity.neighbors(i, neighbors);
        neighbors.push_back(i);

        std::vector<Local_point_2> points(neighbors.size());
        for (std::size_t i = 0; i < neighbors.size(); ++i) {

          CGAL_precondition(neighbors[i] >= 0);
          CGAL_precondition(neighbors[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + neighbors[i]);
          points[i] = m_to_local_converter(get(m_point_map, key));
        }
        CGAL_precondition(points.size() > 0);

        Local_line_2  fitted_line;
        Local_point_2 fitted_centroid;

        // The best fit line will be a line fitted to all region points with its normal being perpendicular to the line.
        #ifndef CGAL_EIGEN2_ENABLED
          m_scores[i] = linear_least_squares_fitting_2(
            points.begin(), points.end(), 
            fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 2>());
        #else 
          m_scores[i] = linear_least_squares_fitting_2(
            points.begin(), points.end(), 
            fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 2>());
        #endif
      }
      // std::cerr << "* done calculating scores\n";
    }

    const Input_range& m_input_range;
    const Point_map m_point_map;
    Input_connectivity& m_connectivity;
    std::vector<std::size_t> m_order;
    std::vector<Local_FT> m_scores;
    const To_local_converter m_to_local_converter;
    // Seed_map m_seed_map(std::vector<std::size_t>());
  };

} // namespace internal
} // namespace Shape_detection


#endif