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

#ifndef CGAL_BARYCENTRIC_COORDINATES_HARMONIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_COORDINATES_HARMONIC_COORDINATES_2_H

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
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>

// [1] Reference: 

namespace CGAL {
namespace Barycentric_coordinates {

  /*! 
    \ingroup PkgBarycentric_coordinates_2WAC

    \brief Harmonic coordinates.

    This class implements 2D harmonic coordinates (  ) and can be used in conjunction
    with `CGAL::Barycentric_coordinates::pointwise_coordinates_2()` to compute
    evaluate harmonic coordinate functions at any point inside a polygon.
    
    Harmonic coordinates are well-defined and non-negative in the closure 
    of any simple polygon.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam Domain
    is a model of `HarmonicCoordinatesDomain`.

    \tparam Solver
    is a model of `HarmonicCoordinatesSolver`.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam VertexMap 
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`.
    
    \cgalModels `PointwiseWeigts_2`
  */
  template<
  typename Polygon,
  typename Domain,
  typename Solver,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Harmonic_coordinates_2 : public Discrete_harmonic_weights_2<Polygon, GeomTraits, VertexMap> {

    public:

      /// \name Types
      /// @{

      /// \cond SKIP_IN_MANUAL 
      using Base = Discrete_harmonic_weights_2<Polygon, GeomTraits, VertexMap>;
      using Polygon_ = Polygon;
      using Domain_ = Domain;
      using Solver_ = Solver;
      using Traits = GeomTraits;
      using Vertex_map = VertexMap;
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

        This class implements the behavior of harmonic coordinates 
        for a 2D query point.

        \param polygon
        An instance of `Polygon` with vertices of a 2D polygon.

        \param domain
        An instance of `Domain` that contains a triangulation of the polygon's
        interior domain.

        \param solver
        An instance of `Solver` that computes the LP factorization 
        of a sparse matrix.

        \param vertex_map
        An instance of `VertexMap` that maps a vertex from `polygon` 
        to `Point_2`.

        \param traits
        An instance of `GeomTraits`.

        \pre `polygon.size() > 3`
      */
      Harmonic_coordinates_2(
        const Polygon& polygon,
        const Domain& domain,
        const Solver& solver,
        const VertexMap vertex_map = VertexMap(),
        const GeomTraits traits = GeomTraits()) :
      m_polygon(polygon),
      m_domain(domain),
      m_solver(solver),
      m_vertex_map(vertex_map),
      m_traits(traits) {
        
        CGAL_precondition(polygon.size() > 3);
      }

      /// @}

      /// \name Access
      /// @{ 

      /*!
        \brief implements `PointwiseWeights_2::operator()()`.
        
        This function fills `coordinates` with harmonic coordinates 
        evaluated at the point `query` with respect to the vertices of the polygon.
        Evaluation is performed by locating a triangle in the `domain` and then
        applying `CGAL::Barycentric_coordinates::triangles_coordinates_2()` to
        harmonic coordinates associated with the vertices of this triangle.

        This function can be called for any 2D point inside the polygon.

        The number of computed weights is equal to the `polygon.size()`.
        
        \tparam OutputIterator
        is an output iterator whose value type is `FT`.

        \param query
        A query point.

        \param weights
        An output iterator that stores the computed coordinates.
      */
      template<typename OutputIterator>
      boost::optional<OutputIterator> operator()(
        const Point_2& query, 
        OutputIterator coordinates) const {

        if (!is_valid_point(query)) {
          internal::get_default(m_polygon.size(), coordinates);
          return;
        }


      }

      /// @}

      /// \name Verifications
      /// @{

      /*!
        \brief implements `PointwiseWeights_2::is_valid_point()`.

        This function checks if a given point `query` is valid to compute harmonic coordinates.
        It returns `true` if and only if `query` is in the polygon's closure.

        \param query
        A query point.

        \return boolean `true` or `false`.
      */
      bool is_valid_point(const Point_2& query) const {
        
        const auto result = internal::locate_wrt_polygon(
          m_polygon, query, 
          m_vertex_map, 
          m_traits);

        const Query_point_location location = (*result).first;
        if (location == Query_point_location::ON_UNBOUNDED_SIDE) 
          return false;
        
        return true;
      }

      /// @}

      /// \name Computation
      /// @{

      /*!
        computes harmonic coordinates at all vertices of the domain.
      */
      void compute() {

      }

      /// @}

      /// \name Memory Management
      /// @{

      /*!
        clears all internal data structures.
      */
      void clear() {

      }

      /*!
        releases all memory that is used internally.
      */
      void release_memory() {

      }

      /// @}

    private:
      
      // Fields.
      const Polygon& m_polygon;
      const Domain& m_domain;
      const Solver& m_solver;
      const Vertex_map m_vertex_map;
      const Traits m_traits;
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_HARMONIC_COORDINATES_2_H
