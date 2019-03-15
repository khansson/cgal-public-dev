// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_COORDINATES_WACHSPRESS_WEIGHTS_2_H
#define CGAL_BARYCENTRIC_COORDINATES_WACHSPRESS_WEIGHTS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL includes.
#include <vector>

// Boost includes.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos. 
// A general construction of barycentric coordinates over convex polygons. 
// Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

namespace CGAL {
namespace Barycentric_coordinates {

  /*! 
    \ingroup PkgBarycentric_coordinates_2WAC

    \brief Generalized Wachspress weights.

    This class implements 2D Wachspress weights ( \cite cgal:bc:fhk-gcbcocp-06, 
    \cite cgal:bc:mlbd-gbcip-02, \cite cgal:bc:w-rfeb-75 ) and can be used in conjunction
    with `CGAL::Barycentric_coordinates::pointwise_coordinates_2()` to compute
    Wachspress coordinates.
    
    Wachspress coordinates are well-defined and non-negative in the closure 
    of a strictly convex polygon.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap 
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`.
    
    \cgalModels `PointwiseWeigts_2`
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Wachspress_weights_2 {

    public:

      /// \name Types
      /// @{

      /// \cond SKIP_IN_MANUAL  
      using Polygon_ = Polygon;
      using Traits = GeomTraits;
      using Vertex_map = VertexMap;

      using Area_2 = typename Traits::Compute_area_2;
      /// \endcond
      
      /// Number type.
      typedef typename GeomTraits::FT FT;

      /// Point type.
      typedef typename GeomTraits::Point_2 Point_2;

      /// @}

      /// \name Initialization
      /// @{
      
      /*!
        \brief initializes all internal data structures.

        This class implements the behavior of Wachspress weights 
        for a 2D query point.

        \param polygon
        An instance of `Polygon` with vertices of a 2D polygon.

        \param algorithm_type
        The type of the algorithm used to compute weights. %Defaults to the max
        precision version.

        \param vertex_map
        An instance of `VertexMap` that maps a vertex from `polygon` 
        to `Point_2`.

        \param traits
        An instance of `GeomTraits`.

        \pre `polygon.size() > 3`
      */
      Wachspress_weights_2(
        const Polygon& polygon,
        const Algorithm_type algorithm_type = Algorithm_type::MAX_PRECISION,
        const VertexMap vertex_map = VertexMap(),
        const GeomTraits traits = GeomTraits()) :
      m_polygon(polygon),
      m_algorithm_type(algorithm_type),
      m_vertex_map(vertex_map),
      m_traits(traits),
      m_area_2(m_traits.compute_area_2_object()) {
        
        CGAL_precondition(polygon.size() > 3);

        A.resize(polygon.size());
        C.resize(polygon.size());
        w.resize(polygon.size());
      }

      /// @}

      /// \name Access
      /// @{ 

      /*!
        \brief implements `PointwiseWeights_2::operator()()`.
        
        This function fills `weights` with generalized Wachspress weights 
        computed at the point `query` with respect to the vertices of the polygon.

        This function can be called for any 2D point.

        The number of computed weights is equal to the `polygon.size()`.
        
        \tparam OutputIterator
        is an output iterator whose value type is `FT`.

        \param query
        A query point.

        \param weights
        An output iterator that stores the computed weights.

        \pre `algorithm_type == MAX_PRECISION || algorithm_type == MAX_SPEED`.
      */
      template<typename OutputIterator>
      boost::optional<OutputIterator> operator()(
        const Point_2& query, 
        OutputIterator weights) const {

        CGAL_precondition(
          m_algorithm_type == Algorithm_type::MAX_PRECISION ||
          m_algorithm_type == Algorithm_type::MAX_SPEED);

        switch(m_algorithm_type) {
          case Algorithm_type::MAX_PRECISION:
            return max_precision_weights(query, weights);
          case Algorithm_type::MAX_SPEED:
            return max_speed_weights(query, weights);
          default:
            return boost::optional<OutputIterator>();
        }
      }

      /// @}

      /// \name Verifications
      /// @{

      /*!
        \brief implements `PointwiseWeights_2::is_valid_point()`.

        This function checks if a given point `query` is valid to compute Wachspress weights.
        It returns `true` if and only if `query` is in the polygon's closure.
        In fact, Wachspress weights are well-defined everywhere in the plane, but normalized
        weights after applying the function `CGAL::Barycentric_coordinates::pointwise_coordinates_2()`
        are not well-defined outside the polygon. Thus this function is used to control the validity
        of the query point location.

        \param query
        A query point.

        \return boolean `true` or `false`.
      */
      bool is_valid_point(const Point_2& query) const {

        const Polygon_type polygon_type = 
        internal::polygon_type(m_polygon, m_vertex_map, m_traits);

        if (polygon_type == Polygon_type::CONCAVE)
          return false;
        
        const auto result = internal::locate_wrt_polygon(
          m_polygon, query, 
          m_vertex_map, 
          m_traits);

        const Query_point_location location = (*result).first;
        if (location == Query_point_location::ON_UNBOUNDED_SIDE) 
          return false;
        
        return true;
      }

      /*!
        \brief implements `PointwiseWeights_2::is_boundary_point()`.
        
        This function checks if a given point `query` is on the polygon's boundary.

        To compute barycentric coordinates for points that belong to the polygon's boundary,
        it is better to use the function `CGAL::Barycentric_coordinates::boundary_coordinates_2()` 
        instead of computing Wachspress weights and then normalizing them. The latter
        case can cause precision problems.

        \param query
        A query point.

        \return an `std::pair`, where the first item in the pair is location
        of the point `query` with respect to the polygon and second item is
        the index of the polygon vertex or edge if `query` belongs to the
        polygon's boundary. It is std::size_t(-1) if it does not.
      */
      boost::optional< std::pair<Query_point_location, std::size_t> > 
      is_boundary_point(const Point_2& query) const {
        return internal::locate_wrt_polygon(
          m_polygon,
          query,
          m_vertex_map,
          m_traits);
      }

      /// @}

    private:
      
      // Fields.
      std::vector<FT> A;
      std::vector<FT> C;
      std::vector<FT> w;

      const Polygon& m_polygon;
      const Algorithm_type m_algorithm_type;
      const Vertex_map m_vertex_map;
      const Traits m_traits;

      const Area_2 m_area_2;

      // Functions.
      template<typename OutputIterator>
      boost::optional<OutputIterator> max_precision_weights(
        const Point_2& query,
        OutputIterator weights) {

        // Get the number of vertices in the polygon.
        const std::size_t n = m_polygon.size();

        // Compute areas A following the area notation from [1]. 
        // Split the loop to make this computation faster.
        A[0] = m_area_2(m_polygon[0], m_polygon[1], query);
        for (std::size_t i = 1; i < n-1; ++i) 
          A[i] = m_area_2m_polygon[i], m_polygon[i+1], query);
        A[n-1] = m_area_2(m_polygon[n-1], m_polygon[0], query);

        // Initialize weights with areas C following the area notation from [1].
        // Then we multiply them by areas A as in the formula (5) from [1]. 
        // We also split the loop.
        w[0] = m_area_2(m_polygon[n-1], m_polygon[0], m_polygon[1]);
        for(std::size_t j = 1; j < n-1; ++j) 
          w[0] *= A[j];

        for(std::size_t i = 1; i < n-1; ++i) {
          w[i] = m_area_2(m_polygon[i-1], m_polygon[i], m_polygon[i+1]);
          
          for (std::size_t j = 0; j < i-1; ++j) w[i] *= A[j];
          for (std::size_t j = i+1; j < n; ++j) w[i] *= A[j];
        }

        w[n-1] = m_area_2(m_polygon[n-2], m_polygon[n-1], m_polygon[0]);
        for (std::size_t j = 0; j < n-2; ++j) 
          w[n-1] *= A[j];

        // Return weights.
        for (std::size_t i = 0; i < n; ++i)
          *(weights++) = w[i]

        return boost::optional<OutputIterator>(weights);
      }

      template<typename OutputIterator>
      boost::optional<OutputIterator> max_speed_weights(
        const Point_2& query,
        OutputIterator weights) {

        // Get the number of vertices in the polygon.
        const std::size_t n = m_polygon.size();

        // Compute areas A and C following the area notation from [1]. 
        // Split the loop to make this computation faster.
        A[0] = m_area_2(m_polygon[0], m_polygon[1], query);
        C[0] = m_area_2(m_polygon[n-1], m_polygon[0], m_polygon[1]);

        for (std::size_t i = 1; i < n-1; ++i) {
          A[i] = m_area_2(m_polygon[i], m_polygon[i+1], query);
          C[i] = m_area_2(m_polygon[i-1], m_polygon[i], m_polygon[i+1]);
        }

        A[n-1] = m_area_2(m_polygon[n-1], m_polygon[0], query);
        C[n-1] = m_area_2(m_polygon[n-2], m_polygon[n-1], m_polygon[0]);

        // Compute unnormalized weights following the formula (28) from [1].
        CGAL_precondition(A[n-1] != FT(0) && A[0] != FT(0));
        *(weights++) = C[0] / (A[n-1] * A[0]);

        for (std::size_t i = 1; i < n-1; ++i) {
          CGAL_precondition(A[i-1] != FT(0) && A[i] != FT(0));
          *(weights++) = C[i] / (A[i-1] * A[i]);
        }

        CGAL_precondition(A[n-2] != FT(0) && A[n-1] != FT(0));
        *(weights++) = C[n-1] / (A[n-2] * A[n-1]);

        // Return weights.
        return boost::optional<OutputIterator>(weights);
      }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_WACHSPRESS_WEIGHTS_2_H
