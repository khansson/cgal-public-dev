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

#ifndef CGAL_BARYCENTRIC_MEAN_VALUE_WEIGHTS_2_H
#define CGAL_BARYCENTRIC_MEAN_VALUE_WEIGHTS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL includes.
#include <vector>
#include <utility>
#include <iterator>

// Boost includes.
#include <boost/optional/optional.hpp>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// [1] Reference: "K. Hormann and M. Floater. 
// Mean value coordinates for arbitrary planar polygons. 
// ACM Transactions on Graphics, 25(4):1424-1441, 2006.".
// [2] Reference: "M. S. Floater, 
// Wachspress and mean value coordinates, 
// to appear in the Proceedings of the 14th International Conference on Approximation Theory, 
// G. Fasshauer and L. L. Schumaker (eds.)."

namespace CGAL {
namespace Barycentric_coordinates {

  /*! 
    \ingroup PkgBarycentric_coordinates_2MVC

    \brief Mean value weights.

    This class implements 2D mean value weights ( \cite cgal:bc:hf-mvcapp-06, 
    \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ) and can be used in conjunction
    with `Barycentric_coordinates::analytic_coordinates_2()` to compute
    mean value coordinates.
    
    Mean value coordinates are well-defined everywhere in the plane and are 
    non-negative in the kernel of a star-shaped polygon.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap 
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \cgalModels `AnalyticWeights_2`
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Mean_value_weights_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Polygon_ = Polygon;
    using GeomTraits_ = GeomTraits;
    using VertexMap_ = VertexMap;

    using Vector_2 = typename GeomTraits::Vector_2;
    using Area_2 = typename GeomTraits::Compute_area_2;
    using Squared_length_2 = typename GeomTraits::Compute_squared_length_2;
    using Scalar_product_2 = typename GeomTraits::Compute_scalar_product_2;
    using Get_sqrt = internal::Get_sqrt<GeomTraits>;
    using Sqrt = typename Get_sqrt::Sqrt;
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

      This class implements the behavior of mean value weights 
      for 2D query points.

      \param polygon
      An instance of `Polygon` with vertices of a simple polygon.

      \param computation_policy
      One of the `Barycentric_coordinates::Computation_policy`.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon` 
      to `Point_2`.

      \param traits
      An instance of `GeomTraits`.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
    */
    Mean_value_weights_2(
      const Polygon& input_polygon,
      const Computation_policy computation_policy 
        = Computation_policy::PRECISE_COMPUTATION_WITH_EDGE_CASES,
      const VertexMap vertex_map = VertexMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_polygon(input_polygon),
    m_computation_policy(computation_policy),
    m_vertex_map(vertex_map),
    m_traits(traits),
    m_area_2(m_traits.compute_area_2_object()),
    m_squared_length_2(m_traits.compute_squared_length_2_object()),
    m_scalar_product_2(m_traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(m_traits))  {
        
      m_polygon.clear();
      m_polygon.reserve(m_input_polygon.size());
      for (const auto& item : m_input_polygon)
        m_polygon.push_back(get(m_vertex_map, item));

      CGAL_precondition(m_polygon.size() >= 3);

      s.resize(m_polygon.size());
      r.resize(m_polygon.size());
      A.resize(m_polygon.size());
      B.resize(m_polygon.size());
      D.resize(m_polygon.size());
      P.resize(m_polygon.size());
      t.resize(m_polygon.size());
      w.resize(m_polygon.size());

      CGAL_precondition(
        CGAL::is_simple_2(m_polygon.begin(), m_polygon.end(), m_traits));
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief implements `AnalyticWeights_2::operator()()`.
        
      This function fills `weights` with mean value weights 
      computed at the `query` point with respect to the vertices of the polygon.
      If `query` belongs to the polygon's boundary, the returned weights are normalized.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param weights
      An output iterator that stores the computed weights.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> operator()(
      const Polygon&,
      const Point_2& query, 
      OutputIterator weights,
      GeomTraits) {

      switch(m_computation_policy) {
        case Computation_policy::PRECISE_COMPUTATION: {
          return max_precision_weights(query, weights);
        }
        case Computation_policy::PRECISE_COMPUTATION_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, weights);
          if (edge_case == internal::Edge_case::BOUNDARY)
            return weights;
          return max_precision_weights(query, weights);
        }
        case Computation_policy::FAST_COMPUTATION: {
          return max_speed_weights(query, weights);
        }
        case Computation_policy::FAST_COMPUTATION_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, weights);
          if (edge_case == internal::Edge_case::BOUNDARY)
            return weights;
          return max_speed_weights(query, weights);
        }
        default: {
          const auto edge_case = verify(query, weights);
          if (edge_case == internal::Edge_case::BOUNDARY)
            return weights;
          return max_precision_weights(query, weights);
        }
      }
      return boost::none;
    }

    /*! 
      This function fills `weights` with mean value weights 
      computed at the `query` point with respect to the vertices of the polygon.
        
      This function calls the function above.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param weights
      An output iterator that stores the computed weights.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> operator()(
      const Point_2& query, 
      OutputIterator weights) {
      
      return operator()(m_input_polygon, query, weights, m_traits);
    }

    /// @}

  private:
      
    // Fields.
    std::vector<Vector_2> s;

    std::vector<FT> r;
    std::vector<FT> A;
    std::vector<FT> B;
    std::vector<FT> D;
    std::vector<FT> P;
    std::vector<FT> t;
    std::vector<FT> w;

    const Polygon& m_input_polygon;
    const Computation_policy m_computation_policy;
    const VertexMap m_vertex_map;
    const GeomTraits m_traits;

    const Area_2 m_area_2;
    const Squared_length_2 m_squared_length_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;
    
    std::vector<Point_2> m_polygon;

    // Functions.
    template<typename OutputIterator>
    internal::Edge_case verify(
      const Point_2& query,
      OutputIterator weights) const {

      const auto result = internal::locate_wrt_polygon_2(
        m_polygon, query, m_traits);

      const Query_point_location location = (*result).first;
      const std::size_t index = (*result).second;

      if (location == Query_point_location::ON_UNBOUNDED_SIDE) 
        return internal::Edge_case::UNBOUNDED;

      if (
        location == Query_point_location::ON_VERTEX ||
        location == Query_point_location::ON_EDGE ) {
        internal::boundary_coordinates_2(
          m_polygon, query, location, index, weights, m_traits);
        return internal::Edge_case::BOUNDARY;
      }

      return internal::Edge_case::INTERIOR;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> max_precision_weights(
      const Point_2& query,
      OutputIterator weights) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute vectors s and its lengths r following the pseudo-code 
      // in the Figure 10 from [1].
      s[0] = m_polygon[0] - query;
      r[0] = m_sqrt(m_squared_length_2(s[0]));

      // Compute areas A and B following the notation from [1] (see Figure 2). 
      // Split the loop to make this computation faster.
      A[0] = m_area_2(m_polygon[0], m_polygon[1], query);
      B[0] = m_area_2(m_polygon[n-1], m_polygon[1], query);

      for (std::size_t i = 1; i < n-1; ++i) {
        s[i] = m_polygon[i] - query;
        r[i] = m_sqrt(m_squared_length_2(s[i]));

        A[i] = m_area_2(m_polygon[i], m_polygon[i+1], query);
        B[i] = m_area_2(m_polygon[i-1], m_polygon[i+1], query);
      }

      s[n-1] = m_polygon[n-1] - query;
      r[n-1] = m_sqrt(m_squared_length_2(s[n-1]));

      A[n-1] = m_area_2(m_polygon[n-1], m_polygon[0], query);
      B[n-1] = m_area_2(m_polygon[n-2], m_polygon[0], query);

      // Following section 4.2 from [2] we denote P_j = r_j*r_{j+1} + dot_product(d_j, d_{j+1}).
      // Vector s_i from [1] corresponds to that one with the name d_i in [2].
      for (std::size_t j = 0; j < n-1; ++j)
        P[j] = CGAL::max(r[j] * r[j+1] + m_scalar_product_2(s[j], s[j+1]), FT(0));
      P[n-1] = CGAL::max(r[n-1] * r[0] + m_scalar_product_2(s[n-1], s[0]), FT(0));

      // Compute mean value weights using the formula (16) from [2].
      // Since the formula (16) always gives positive values, 
      // we have to add a proper sign to all the weight functions.
      w[0] = r[n-1]*r[1] - m_scalar_product_2(s[n-1], s[1]);
      for (std::size_t j = 1; j < n-1; ++j) 
        w[0] *= P[j];
      w[0] = sign_of_weight(A[n-1], A[0], B[0]) * m_sqrt(w[0]);

      for (std::size_t i = 1; i < n-1; ++i) {
        w[i] = r[i-1] * r[i+1] - m_scalar_product_2(s[i-1], s[i+1]);
          
        for (std::size_t j = 0; j < i-1; ++j) 
          w[i] *= P[j];
          
        for(std::size_t j = i+1; j < n; ++j) 
          w[i] *= P[j];
          
        w[i] = sign_of_weight(A[i-1], A[i], B[i]) * m_sqrt(w[i]);
      }

      w[n-1] = r[n-2] * r[0] - m_scalar_product_2(s[n-2], s[0]);
      for (std::size_t j = 0; j < n-2; ++j) 
        w[n-1] *= P[j];
      w[n-1] = sign_of_weight(A[n-2], A[n-1], B[n-1]) * m_sqrt(w[n-1]);

      // Return weights.
      for (std::size_t i = 0; i < n; ++i)
        *(weights++) = w[i];

      return weights;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> max_speed_weights(
      const Point_2& query,
      OutputIterator weights) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute vectors s following the pseudo-code in the Figure 10 from [1].
      for (std::size_t i = 0; i < n; ++i) 
        s[i] = m_polygon[i] - query;

      // Compute lengths r, areas A, and dot products D following the pseudo-code 
      // in the Figure 10 from [1].
      // Split the loop to make this computation faster.
      r[0] = m_sqrt(m_squared_length_2(s[0]));
      A[0] = m_area_2(m_polygon[0], m_polygon[1], query);
      D[0] = m_scalar_product_2(s[0], s[1]);

      for (std::size_t i = 1; i < n-1; ++i) {
        r[i] = m_sqrt(m_squared_length_2(s[i]));
        A[i] = m_area_2(m_polygon[i], m_polygon[i+1], query);
        D[i] = m_scalar_product_2(s[i], s[i+1]);
      }

      r[n-1] = m_sqrt(m_squared_length_2(s[n-1]));
      A[n-1] = m_area_2(m_polygon[n-1], m_polygon[0], query);
      D[n-1] = m_scalar_product_2(s[n-1], s[0]);

      // Compute intermediate values t using the formulas from slide 19 here
      // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf
      for (std::size_t i = 0; i < n-1; ++i) {
        CGAL_precondition((r[i]*r[i+1] + D[i]) != FT(0));
        t[i] = A[i] / (r[i]*r[i+1] + D[i]);
      }

      CGAL_precondition((r[n-1]*r[0] + D[n-1]) != FT(0));
      t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

      // Compute mean value weights using the same pseudo-code as before.
      CGAL_precondition(r[0] != FT(0));
      *(weights++) = (t[n-1] + t[0]) / r[0];

      for (std::size_t i = 1; i < n-1; ++i) {
        CGAL_precondition(r[i] != FT(0));
        *(weights++) = (t[i-1] + t[i]) / r[i];
      }

      CGAL_precondition(r[n-1] != FT(0));
      *(weights++) = (t[n-2] + t[n-1]) / r[n-1];

      // Return weights.
      return weights;
    }

    // Return the sign of a mean value weight function.
    // We can have 3 different values: 0 if the weight = 0, -1 
    // if the weight is negative, and +1 if the weight is positive.
    FT sign_of_weight(const FT& A_prev, const FT& A, const FT& B) const {
        
      if (A_prev > FT(0) && A > FT(0) && B <= FT(0)) return FT(1);
      if (A_prev < FT(0) && A < FT(0) && B >= FT(0)) return FT(-1);
      if (B > FT(0)) return FT(1);
      if (B < FT(0)) return FT(-1);

      return FT(0);
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_MEAN_VALUE_WEIGHTS_2_H
