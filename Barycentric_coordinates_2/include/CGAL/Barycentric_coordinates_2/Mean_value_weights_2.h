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

#ifndef CGAL_BARYCENTRIC_COORDINATES_MEAN_VALUE_WEIGHTS_2_H
#define CGAL_BARYCENTRIC_COORDINATES_MEAN_VALUE_WEIGHTS_2_H

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
    \ingroup PkgBarycentric_coordinates_2WAC

    \brief Generalized mean value weights.

    This class implements 2D mean value weights ( \cite cgal:bc:fhk-gcbcocp-06, 
    \cite cgal:pp-cdmsc-93, \cite cgal:bc:eddhls-maam-95 ) and can be used in conjunction
    with `CGAL::Barycentric_coordinates::pointwise_coordinates_2()` to compute
    mean value coordinates.
    
    Mean value coordinates are well-defined everywhere in the plane and are 
    non-negative in the kernel of a star-shaped polygon.

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
  class Mean_value_weights_2 {

    public:

      /// \name Types
      /// @{

      /// \cond SKIP_IN_MANUAL  
      using Polygon_ = Polygon;
      using Traits = GeomTraits;
      using Vertex_map = VertexMap;

      using Area_2 = typename Traits::Compute_area_2;
      using Squared_length_2 = typename Traits::Compute_squared_length_2;
      using Scalar_product_2 = typename Traits::Compute_scalar_product_2;
      using Get_sqrt = internal::Get_sqrt<Traits>;
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
      Mean_value_weights_2(
        const Polygon& polygon,
        const Algorithm_type algorithm_type = Algorithm_type::MAX_PRECISION,
        const VertexMap vertex_map = VertexMap(),
        const GeomTraits traits = GeomTraits()) :
      m_polygon(polygon),
      m_algorithm_type(algorithm_type),
      m_vertex_map(vertex_map),
      m_traits(traits),
      m_area_2(m_traits.compute_area_2_object()),
      m_squared_length_2(m_traits.compute_squared_length_2_object()),
      m_scalar_product_2(m_traits.compute_scalar_product_2_object()),
      m_sqrt(Get_sqrt::sqrt_object(m_traits)) {
        
        CGAL_precondition(polygon.size() > 3);

        s.resize(polygon.size());
        r.resize(polygon.size());
        A.resize(polygon.size());
        B.resize(polygon.size());
        D.resize(polygon.size());
        P.resize(polygon.size());
        t.resize(polygon.size());
        w.resize(polygon.size());
      }

      /// @}

      /// \name Access
      /// @{ 

      /*!
        \brief implements `PointwiseWeights_2::operator()()`.
        
        This function fills `weights` with generalized mean value weights 
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

        This function checks if a given point `query` is valid to compute mean value weights.
        It always returns `true`.

        \param query
        A query point.

        \return boolean `true` or `false`.
      */
      bool is_valid_point(const Point_2& query) const {
        return true;
      }

      /*!
        \brief implements `PointwiseWeights_2::is_boundary_point()`.
        
        This function checks if a given point `query` is on the polygon's boundary.

        To compute barycentric coordinates for points that belong to the polygon's boundary,
        it is better to use the function `CGAL::Barycentric_coordinates::boundary_coordinates_2()` 
        instead of computing mean value weights and then normalizing them. The latter
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
      std::vector<Point_2> s;

      std::vector<FT> r;
      std::vector<FT> A;
      std::vector<FT> B;
      std::vector<FT> D;
      std::vector<FT> P;
      std::vector<FT> t;
      std::vector<FT> w;

      const Polygon& m_polygon;
      const Algorithm_type m_algorithm_type;
      const Vertex_map m_vertex_map;
      const Traits m_traits;

      const Area_2 m_area_2;
      const Squared_length_2 m_squared_length_2;
      const Scalar_product_2 m_scalar_product_2;
      const Sqrt m_sqrt;

      // Functions.
      template<typename OutputIterator>
      boost::optional<OutputIterator> max_precision_weights(
        const Point_2& query,
        OutputIterator weights) {

      }

      template<typename OutputIterator>
      boost::optional<OutputIterator> max_speed_weights(
        const Point_2& query,
        OutputIterator weights) {

      }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_MEAN_VALUE_WEIGHTS_2_H
