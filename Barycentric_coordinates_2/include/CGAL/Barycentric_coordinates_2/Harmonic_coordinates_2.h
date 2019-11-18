// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_BARYCENTRIC_HARMONIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_HARMONIC_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL includes.
#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <iterator>

// Boost includes.
#include <boost/optional/optional.hpp>

// Eigen includes.
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_weights_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// [1] Reference: 

namespace CGAL {
namespace Barycentric_coordinates {

  /*! 
    \ingroup PkgBarycentric_coordinates_2HAC

    \brief Harmonic coordinates.

    This class implements 2D harmonic coordinates and can be used in conjunction
    with `Barycentric_coordinates::analytic_coordinates_2()` to evaluate 
    harmonic coordinate functions at any point inside a polygon.
    
    Harmonic coordinates are well-defined and non-negative in the closure 
    of any simple polygon.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam Triangulation
    is a model of `HarmonicCoordinatesDomain_2`.

    \tparam Solver
    is a model of `HarmonicCoordinatesSolver_2`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap 
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \cgalModels `AnalyticWeights_2`
  */
  template<
  typename Polygon,
  typename Triangulation,
  typename Solver,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Harmonic_coordinates_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL 
    using Polygon_ = Polygon;
    using Triangulation_ = Triangulation;
    using Solver_ = Solver;
    using GeomTraits_ = GeomTraits;
    using VertexMap_ = VertexMap;
    using Face_handle = typename Triangulation::Face_handle;
    /// \endcond
      
    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// Vertex handle type.
    typedef typename Triangulation::Vertex_handle Vertex_handle;

    /// \cond SKIP_IN_MANUAL 
    using MatrixFT  = Eigen::SparseMatrix<FT>;
    using VectorFT  = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
    using TripletFT = Eigen::Triplet<FT>;

    using Discrete_harmonic = 
      CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2<
        Polygon, GeomTraits, VertexMap>;
    using Policy = CGAL::Barycentric_coordinates::Computation_policy;
    using Location = std::pair<Query_point_location, std::size_t>;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{
      
    /*!
      \brief initializes all internal data structures.

      This class implements the behavior of harmonic coordinates 
      for 2D query points.

      \param input_polygon
      An instance of `Polygon` with vertices of a simple polygon.

      \param triangulation
      An instance of `Triangulation`.

      \param solver
      An instance of `Solver`.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon` 
      to `Point_2`.

      \param traits
      An instance of `GeomTraits`.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
    */
    Harmonic_coordinates_2(
      const Polygon& input_polygon,
      const Triangulation& triangulation,
      const Solver& solver,
      const VertexMap vertex_map = VertexMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_polygon(input_polygon),
    m_triangulation(triangulation),
    m_solver(solver),
    m_vertex_map(vertex_map),
    m_traits(traits) { 

      m_polygon.clear();
      m_polygon.reserve(m_input_polygon.size());
      for (const auto& item : m_input_polygon)
        m_polygon.push_back(get(m_vertex_map, item));

      CGAL_precondition(
        m_polygon.size() >= 3);
      CGAL_precondition(
        CGAL::is_simple_2(m_polygon.begin(), m_polygon.end(), m_traits));
      CGAL_precondition(
        m_triangulation.number_of_vertices() >= m_polygon.size());
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief implements `AnalyticWeights_2::operator()()`.
        
      This function fills `coordinates` with harmonic coordinates 
      evaluated at the point `query` with respect to the vertices of the polygon.
      Evaluation is performed by locating a triangle in the `triangulation` and then
      interpolating harmonic coordinates within this triangle.
        
      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param coordinates
      An output iterator that stores the computed coordinates.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> operator()(
      const Polygon&,
      const Point_2& query, 
      OutputIterator coordinates,
      GeomTraits) {

      const auto fh = m_triangulation.locate(query);
      if (!fh->is_in_domain()) {
        internal::get_default(m_polygon.size(), coordinates);
        return coordinates;
      }

      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();

      std::vector<FT> b;
      b.reserve(3);

      CGAL::Barycentric_coordinates::triangle_coordinates_2(
        p0, p1, p2, query, std::back_inserter(b), m_traits);

      CGAL_assertion(b.size() == 3);
      CGAL_assertion(m_coordinates.find(fh->vertex(0)) != m_coordinates.end());
      CGAL_assertion(m_coordinates.find(fh->vertex(1)) != m_coordinates.end());
      CGAL_assertion(m_coordinates.find(fh->vertex(2)) != m_coordinates.end());

      const auto& hm0 = m_coordinates.at(fh->vertex(0));
      const auto& hm1 = m_coordinates.at(fh->vertex(1));
      const auto& hm2 = m_coordinates.at(fh->vertex(2));

      std::vector<FT> result(m_polygon.size(), FT(0));
      for (std::size_t i = 0; i < m_polygon.size(); ++i)
        result[i] = hm0[i] * b[0] + hm1[i] * b[1] + hm2[i] * b[2];
      
      for (std::size_t i = 0; i < result.size(); ++i)
        *(coordinates++) = result[i];
      
      return coordinates;
    }

    /*! 
      This function fills `coordinates` with harmonic coordinates 
      evaluated at the point `query` with respect to the vertices of the polygon.
      Evaluation is performed by locating a triangle in the `triangulation` and then
      interpolating harmonic coordinates within this triangle.
        
      This function calls the function above.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param coordinates
      An output iterator that stores the computed coordinates.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> operator()(
      const Point_2& query, 
      OutputIterator coordinates) {
      
      return operator()(m_input_polygon, query, coordinates, m_traits);
    }

    /*!
      \brief returns harmonic coordinates computed at the given triangulation vertex.
        
      This function fills `coordinates` with harmonic coordinates 
      computed at the vertex `vh` of the input triangulation.
        
      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param vh
      A vertex handle.

      \param coordinates
      An output iterator that stores the computed coordinates.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> coordinates(
      const Vertex_handle vh, 
      OutputIterator coordinates) const {

      const auto& bs = m_coordinates.at(vh);
      for (const FT& b : bs)
        *(coordinates++) = b;
      return coordinates;
    }

    /*!
      \brief returns harmonic coordinates computed at the triangulation vertices.
        
      This function fills `coordinates` with harmonic coordinates 
      computed at all the vertices of the input triangulation.
        
      \tparam OutputIterator
      is an output iterator whose value type is `std::vector<FT>`.

      \param coordinates
      An output iterator that stores the computed coordinates.
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> coordinates(
      OutputIterator coordinates) const {

      for (const auto& item : m_coordinates) {
        const auto& b = item.second;
        *(coordinates++) = b;
      }
    }

    /// @}

    /// \name Computation
    /// @{

    /*!
      computes harmonic coordinates at all the vertices of the input triangulation.
    */
    void compute() {

      const std::size_t n = m_polygon.size();
      const std::size_t N = m_triangulation.number_of_vertices();

      std::vector<Location> tags;
      tags.reserve(N);

      for (auto vh = m_triangulation.finite_vertices_begin();
      vh != m_triangulation.finite_vertices_end(); ++vh) {
        const auto& query = vh->point();
        tags.push_back(*(internal::locate_wrt_polygon_2(
          m_polygon, query, m_traits)));
      }
      CGAL_assertion(tags.size() == N);

      std::vector<std::size_t> indices;
      indices.reserve(N);

      std::size_t out_count = 0;
      std::size_t numB = 0, numI = 0;

      for (auto vh = m_triangulation.finite_vertices_begin();
      vh != m_triangulation.finite_vertices_end(); ++vh, ++out_count) {
        if (
          tags[out_count].first == Query_point_location::ON_VERTEX ||
          tags[out_count].first == Query_point_location::ON_EDGE ) 
          indices.push_back(numB++);
        else if (
          tags[out_count].first == Query_point_location::ON_BOUNDED_SIDE)
          indices.push_back(numI++);
      }
      CGAL_assertion(indices.size() <= N);

      VectorFT boundary(numB, n);
      MatrixFT A(numI, numI);

      VectorFT x = VectorFT::Zero(numI, n);
      VectorFT b = VectorFT::Zero(numI, n);

      out_count = 0;
      std::size_t in_count = 0;
      std::vector<FT> lambda; 
      lambda.reserve(n);

      std::vector<FT> empty_vec;
      internal::get_default(n, std::back_inserter(empty_vec));

      for (auto vh = m_triangulation.finite_vertices_begin();
      vh != m_triangulation.finite_vertices_end(); ++vh, ++out_count) {
        const auto& query = vh->point();
        m_coordinates[vh] = empty_vec;
        
        if (
          tags[out_count].first == Query_point_location::ON_VERTEX ||
          tags[out_count].first == Query_point_location::ON_EDGE ) {
          
          lambda.clear();
          internal::boundary_coordinates_2(
            m_polygon, query, tags[out_count].first, tags[out_count].second,
            std::back_inserter(lambda), m_traits);
          
          for (std::size_t j = 0; j < n; ++j) 
            boundary(indices[in_count], j) = lambda[j];
          ++in_count;
        } else if (
          tags[out_count].first == Query_point_location::ON_BOUNDED_SIDE)
          ++in_count;
      }

      std::vector<TripletFT> tripletList;
      tripletList.reserve(numI * 7);

      out_count = 0;
      for (auto vh = m_triangulation.finite_vertices_begin();
      vh != m_triangulation.finite_vertices_end(); ++vh, ++out_count) {
        const auto& query = vh->point();
        
        if (
          tags[out_count].first == Query_point_location::ON_BOUNDED_SIDE) {

          // std::vector<int> neighbours;
          // _mesh.getRing(i, neighbours);

          // const std::size_t nn = neighbours.size();
          // std::vector<double> alphaCot(nn), betaCot(nn);

          // for (std::size_t j = 0; j < nn; ++j) {
          //   const size_t jp = (j + 1) % nn;

          //   VertexR2 s1 = p[i] - p[neighbours[j]];
          //   VertexR2 s2 = p[neighbours[jp]] - p[neighbours[j]];

          //   alphaCot[j] = cotangent(s2, s1);

          //   s1 = p[neighbours[j]] - p[neighbours[jp]];
          //   s2 = p[i] - p[neighbours[jp]];

          //   betaCot[j] = cotangent(s2, s1);
          // }

          // double W = 0.0;
          // for (std::size_t j = 0; j < nn; ++j) {

          //   const std::size_t jp  = (j + 1) % nn;
          //   const std::size_t idx = neighbours[jp];

          //   const double w = -(alphaCot[j] + betaCot[jp]);
          //   W -= w;

          //   if (p[idx].type != INTERIOR) {
          //     for (std::size_t k = 0; k < n; ++k)
          //       b(indices[i], k) -= boundary(indices[idx], k) * w;
          //   } else {
          //     tripletList.push_back(T(indices[i], indices[idx], w));
          //   }
          // }
          // tripletList.push_back(T(indices[i], indices[i], W));
        }
      }

      A.setFromTriplets(
        tripletList.begin(), tripletList.end());
      A.makeCompressed();

      /* solveLinearSystem(A, b, x); */

      for (std::size_t i = 0; i < n; ++i) {        
        out_count = 0; in_count = 0;

        for (auto vh = m_triangulation.finite_vertices_begin();
        vh != m_triangulation.finite_vertices_end(); ++vh, ++out_count) {
          if (
            tags[out_count].first == Query_point_location::ON_VERTEX ||
            tags[out_count].first == Query_point_location::ON_EDGE ) {
                    
            auto& vec = m_coordinates.at(vh);
            vec[i] = boundary(indices[in_count], i);
            ++in_count;

          } else if (
            tags[out_count].first == Query_point_location::ON_BOUNDED_SIDE) {
                        
            auto& vec = m_coordinates.at(vh);
            vec[i] = x(indices[in_count], i);
            ++in_count;
          }
        }
      }
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      clears all internal data structures.
    */
    void clear() {
      m_coordinates.clear();
    }

    /*!
      releases all memory that is used internally.
    */
    void release_memory() {
      m_coordinates.shrink_to_fit();
    }

    /// @}

  private:
      
    // Fields.
    const Polygon& m_input_polygon;
    const Triangulation& m_triangulation;
    const Solver& m_solver;
    const VertexMap m_vertex_map;
    const GeomTraits m_traits;

    std::vector<Point_2> m_polygon;
    std::map<Vertex_handle, std::vector<FT> > m_coordinates;

    // Function that solves the linear system.
    void solve_linear_system(
      const MatrixFT& A, const VectorFT& b, VectorFT& x) const {

      m_solver.compute(A);
      x = m_solver.solve(b);
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_HARMONIC_COORDINATES_2_H
