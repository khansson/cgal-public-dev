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
#include <map>
#include <vector>
#include <utility>
#include <iterator>

// Boost includes.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/pointwise_coordinates_2.h>

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

    \tparam Triangulation
    is any type that inherits from `CGAL::Triangulation_2`.

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
  typename Triangulation,
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
    using Triangulation_ = Triangulation;
    using Solver_ = Solver;
    using Traits = GeomTraits;
    using Vertex_map = VertexMap;

    using Vertex_handle = typename Triangulation::Vertex_handle;
    using Face_handle = typename Triangulation::Face_handle;
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

      \param triangulation
      An instance of `Triangulation` that contains a triangulation of the polygon's
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
      const Triangulation& triangulation,
      const Solver& solver,
      const VertexMap vertex_map = VertexMap(),
      const GeomTraits traits = GeomTraits()) :
    m_polygon(polygon),
    m_triangulation(triangulation),
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
      Evaluation is performed by locating a triangle in the `triangulation` and then
      applying `CGAL::Barycentric_coordinates::triangles_coordinates_2()` to
      harmonic coordinates associated with the vertices of this triangle.

      This function can be called for any 2D point inside the polygon.

      The number of computed weights is equal to the `polygon.size()`.
        
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
      OutputIterator coordinates) const {

      if (!is_valid_point(query)) {
        internal::get_default(m_polygon.size(), coordinates);
        return;
      }

      const auto fh = m_triangulation.locate(query);
      const Point_2& p1 = fh->vertex(0)->point();
      const Point_2& p2 = fh->vertex(1)->point();
      const Point_2& p3 = fh->vertex(2)->point();

      std::vector<FT> b;
      b.reserve(3);

      CGAL::triangle_coordinates_2(p1, p2, p3, query, 
      std::back_inserter(b), m_traits);
      CGAL_assertion(b.size() == 3);

      const auto hm0 = m_coordinates.at(fh->vertex(0));
      const auto hm1 = m_coordinates.at(fh->vertex(1));
      const auto hm2 = m_coordinates.at(fh->vertex(2));

      std::vector<FT> result(m_polygon.size(), FT(0));
      for (std::size_t i = 0; i < m_polygon.size(); ++i)
        result[i] = hm0[i] * b[0] + hm1[i] * b[1] + hm2[i] * b[2];
      
      for (std::size_t i = 0; i < result.size(); ++i)
        *(coordinates++) = result[i];
      
      return boost::optional<OutputIterator>(coordinates);
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
      computes harmonic coordinates at all vertices of the triangulation.
    */
    void compute() {

      // Boundary and interior.
      const std::size_t n = _v.size();
      const std::size_t N = _mesh.numVertices();

      CGAL_assertion(N != 0);
      bb.resize(N);

      std::vector<VertexR2> &p = _mesh.vertices();
      std::vector<std::size_t> indices(N, 0);

      std::size_t numB = 0, numI = 0;
      for (std::size_t i = 0; i < N; ++i) {

        if (p[i].type != INTERIOR) 
          indices[i] = numB++;
        else 
          indices[i] = numI++;
      }

      MatrixXd boundary(numB, n);
      Eigen::SparseMatrix<double> A(numI, numI);

      MatrixXd x = MatrixXd::Zero(numI, n);
      MatrixXd b = MatrixXd::Zero(numI, n);

      for (std::size_t i = 0; i < N; ++i) {

        p[i].b().clear();
        p[i].b().resize(n, 0.0);

        bb[i].clear();
        bb[i].resize(n, 0.0);

        if (p[i].type != INTERIOR) {
          computeBoundaryCoordinates(p[i]);
          for (std::size_t j = 0; j < n; ++j) 
            boundary(indices[i], j) = p[i].b()[j];
        }
      }

      typedef Eigen::Triplet<double> T;
      std::vector<T> tripletList;
      tripletList.reserve(numI * 7);

      for (std::size_t i = 0; i < N; ++i) {
        if (p[i].type == INTERIOR) {

          std::vector<int> neighbours;
          _mesh.getRing(i, neighbours);

          const std::size_t nn = neighbours.size();
          std::vector<double> alphaCot(nn), betaCot(nn);

          for (std::size_t j = 0; j < nn; ++j) {
            const size_t jp = (j + 1) % nn;

            VertexR2 s1 = p[i] - p[neighbours[j]];
            VertexR2 s2 = p[neighbours[jp]] - p[neighbours[j]];

            alphaCot[j] = cotangent(s2, s1);

            s1 = p[neighbours[j]] - p[neighbours[jp]];
            s2 = p[i] - p[neighbours[jp]];

            betaCot[j] = cotangent(s2, s1);
          }

          double W = 0.0;
          for (std::size_t j = 0; j < nn; ++j) {

            const std::size_t jp  = (j + 1) % nn;
            const std::size_t idx = neighbours[jp];

            const double w = -(alphaCot[j] + betaCot[jp]);
            W -= w;

            if (p[idx].type != INTERIOR) {
              for (std::size_t k = 0; k < n; ++k)
                b(indices[i], k) -= boundary(indices[idx], k) * w;
            } else {
              tripletList.push_back(T(indices[i], indices[idx], w));
            }
          }
          tripletList.push_back(T(indices[i], indices[i], W));
        }
      }

      A.setFromTriplets(tripletList.begin(), tripletList.end());
      A.makeCompressed();
      solveLinearSystem(A, b, x);

      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < N; ++j) {            
          if (p[j].type != INTERIOR) {
                    
            p[j].b()[i] = boundary(indices[j], i);
            bb[j][i] = p[j].b()[i];

          } else {
                        
            p[j].b()[i] = x(indices[j], i);
            bb[j][i] = p[j].b()[i];
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
    const Polygon& m_polygon;
    const Triangulation& m_triangulation;
    const Solver& m_solver;
    const Vertex_map m_vertex_map;
    const Traits m_traits;

    std::map<Vertex_handle, std::vector<FT> > m_coordinates;

    // Function that solves the linear system.
    void solve_linear_system(
      const Eigen::SparseMatrix<double>& A, 
      const MatrixXd& b, 
      MatrixXd& x) const {

      // LDLT solver.
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldlt;

      ldlt.compute(A);
      x = ldlt.solve(b);
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_HARMONIC_COORDINATES_2_H
