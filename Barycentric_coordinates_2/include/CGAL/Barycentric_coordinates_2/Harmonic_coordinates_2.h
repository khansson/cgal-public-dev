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

// Eigen includes.
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// [1] Reference: "P. Joshi, M. Meyer, T. DeRose, B. Green, and T. Sanocki. 
// Harmonic coordinates for character articulation.
// ACM Transactions on Graphics, 26(3):71:1-9, 2007.".

namespace CGAL {
namespace Barycentric_coordinates {

  /*! 
    \ingroup PkgBarycentric_coordinates_2Classes

    \brief Harmonic coordinates.

    This class implements 2D harmonic coordinates and can be used in conjunction
    with `Barycentric_coordinates::analytic_coordinates_2()` to evaluate 
    harmonic coordinate functions at any point inside a polygon.
    
    Harmonic coordinates are well-defined and non-negative in the closure 
    of any simple polygon.

    Internally, the `Barycentric_coordinates::Discrete_harmonic_weights_2`
    are used.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam Domain
    is a model of `DiscretizedDomain_2`. For the moment, we only support domains
    that have elements, which are triangles.

    \tparam Solver
    is a model of `SparseLinearSolver_2`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap 
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \cgalModels `AnalyticWeights_2`
  */
  template<
  typename Polygon,
  typename Domain,
  typename Solver,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Harmonic_coordinates_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL 
    using Polygon_    = Polygon;
    using Domain_     = Domain;
    using Solver_     = Solver;
    using GeomTraits_ = GeomTraits;
    using VertexMap_  = VertexMap;
    /// \endcond
      
    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// \cond SKIP_IN_MANUAL 
    using Vector_2 = typename GeomTraits::Vector_2;

    using MatrixFT  = Eigen::SparseMatrix<FT>;
    using VectorFT  = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
    using TripletFT = Eigen::Triplet<FT>;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{
      
    /*!
      \brief initializes all internal data structures.

      This class implements the behavior of harmonic coordinates 
      for 2D query points.

      \param polygon
      An instance of `Polygon` with vertices of a simple polygon.

      \param domain
      An instance of `Domain`.

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
      const Polygon& polygon,
      const Domain& domain,
      Solver& solver,
      const VertexMap vertex_map = VertexMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_polygon(polygon),
    m_domain(domain),
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
        m_domain.number_of_vertices() >= m_polygon.size());

      m_b.reserve(3);
      m_element.reserve(3);
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief implements `AnalyticWeights_2::operator()()`.
        
      This function fills `coordinates` with harmonic coordinates 
      evaluated at the point `query` with respect to the vertices of the polygon.
      Evaluation is performed by locating an element in the `domain` that contains
      `query` and then interpolating harmonic coordinates within this element.
        
      If query is not inside `domain`, all coordinates are set to zero. If the
      located element has more than 3 vertices, all coordinates are set to zero.

      Internally, `Barycentric_coordinates::triangle_coordinates_2` are used.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param coordinates
      An output iterator that stores the computed coordinates.

      \pre `element.size() == 3`

      \warning `compute()` should be called before calling this method!
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> operator()(
      const Polygon&,
      const Point_2& query, 
      OutputIterator coordinates,
      GeomTraits) {

      const std::size_t n = m_polygon.size();

      m_element.clear();
      const bool is_in_domain = m_domain.locate(query, m_element);
      if (!is_in_domain || m_element.size() > 3) {
        internal::get_default(n, coordinates);
        return coordinates;
      }

      const std::size_t i0 = m_element[0];
      const std::size_t i1 = m_element[1];
      const std::size_t i2 = m_element[2];

      CGAL_assertion(i0 >= 0 && i0 < m_domain.number_of_vertices());
      CGAL_assertion(i1 >= 0 && i1 < m_domain.number_of_vertices());
      CGAL_assertion(i2 >= 0 && i2 < m_domain.number_of_vertices());

      const auto& p0 = m_domain.vertex(i0);
      const auto& p1 = m_domain.vertex(i1);
      const auto& p2 = m_domain.vertex(i2);

      m_b.clear();
      CGAL::Barycentric_coordinates::internal::planar_coordinates_2(
        p0, p1, p2, query, std::back_inserter(m_b), m_traits);
      CGAL_assertion(m_b.size() == 3);
      
      CGAL_assertion(
        m_coordinates.size() == m_domain.number_of_vertices());
      const auto& hm0 = m_coordinates.at(i0);
      const auto& hm1 = m_coordinates.at(i1);
      const auto& hm2 = m_coordinates.at(i2);

      std::vector<FT> result(n, FT(0));
      for (std::size_t k = 0; k < n; ++k)
        *(coordinates++) = hm0[k] * m_b[0] + hm1[k] * m_b[1] + hm2[k] * m_b[2];
      return coordinates;
    }

    /*! 
      This function fills `coordinates` with harmonic coordinates 
      evaluated at the point `query` with respect to the vertices of the polygon.
      Evaluation is performed by locating an element in the `domain` that contains
      `query` and then interpolating harmonic coordinates within this element.
        
      This function calls the generic function above.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param coordinates
      An output iterator that stores the computed coordinates.

      \warning `compute()` should be called before calling this method!
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> operator()(
      const Point_2& query, 
      OutputIterator coordinates) {
      
      return operator()(m_input_polygon, query, coordinates, m_traits);
    }

    /*!
      \brief fills `coordinates` with harmonic coordinates computed at the 
      vertex of the input domain with the index `query`.
        
      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A vertex index.

      \param coordinates
      An output iterator that stores the computed coordinates.

      \pre `query >= 0 && query < domain.number_of_vertices()`

      \warning `compute()` should be called before calling this method!
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> coordinates(
      const std::size_t query, 
      OutputIterator coordinates) const {

      CGAL_precondition(
        query >= 0 && query < m_domain.number_of_vertices());
      CGAL_precondition(
        m_coordinates.size() == m_domain.number_of_vertices());

      const auto& bs = m_coordinates.at(query);
      for (const FT& b : bs)
        *(coordinates++) = b;
      return coordinates;
    }

    /*!
      \brief fills `coordinates` with harmonic coordinates computed at all the 
      vertices of the input domain.
        
      \tparam OutputIterator
      is an output iterator whose value type is `std::vector<FT>`.

      \param coordinates
      An output iterator that stores the computed coordinates.

      \warning `compute()` should be called before calling this method!
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> coordinates(
      OutputIterator coordinates) const {

      CGAL_assertion(
        m_coordinates.size() == m_domain.number_of_vertices());
      for (const auto& b : m_coordinates)
        *(coordinates++) = b;
    }

    /// @}

    /// \name Computation
    /// @{

    /*!
      computes harmonic coordinates at all the vertices of the input domain.
    */
    void compute() {

      const std::size_t n = m_polygon.size();
      const std::size_t N = m_domain.number_of_vertices();

      std::vector<std::size_t> indices;
      indices.reserve(N);

      std::size_t numB = 0, numI = 0;
      for (std::size_t i = 0; i < N; ++i) {
        if (m_domain.is_on_boundary(i)) {
          indices.push_back(numB); ++numB;
        } else {
          indices.push_back(numI); ++numI;
        }
      }

      VectorFT boundary(numB, n);
      MatrixFT A(numI, numI);
      VectorFT x = VectorFT::Zero(numI, n);
      VectorFT b = VectorFT::Zero(numI, n);

      std::vector<FT> empty_vec;
      empty_vec.reserve(n);
      internal::get_default(
        n, std::back_inserter(empty_vec));

      std::vector<FT> lambda; 
      lambda.reserve(n);

      m_coordinates.clear();
      m_coordinates.reserve(N);

      for (std::size_t i = 0; i < N; ++i) {
        m_coordinates.push_back(empty_vec);

        if (m_domain.is_on_boundary(i)) {
          const auto& query = m_domain.vertex(i);
          const auto edge_found = internal::get_edge_index(
            m_polygon, query, m_traits);
          assert(edge_found);

          const auto location = (*edge_found).first;
          const auto index = (*edge_found).second;
          
          lambda.clear();
          internal::boundary_coordinates_2(
            m_polygon, query, location, index, 
            std::back_inserter(lambda), m_traits);

          for (std::size_t k = 0; k < n; ++k)
            boundary(indices[i], k) = lambda[k];
        }
      }

      std::vector<TripletFT> triplet_list;
      triplet_list.reserve(numI * 7);

      std::vector<std::size_t> neighbors;
      std::vector<FT> alpha_cot, beta_cot;

      for (std::size_t i = 0; i < N; ++i) {
        if (!m_domain.is_on_boundary(i)) {
          const auto& query = m_domain.vertex(i);

          neighbors.clear();
          m_domain(i, neighbors);
          const std::size_t nn = neighbors.size();
          
          alpha_cot.clear(); beta_cot.clear();
          for (std::size_t j = 0; j < nn; ++j) {
            const std::size_t jp = (j + 1) % nn;

            const auto& p1 = m_domain.vertex(neighbors[j]);
            const auto& p2 = m_domain.vertex(neighbors[jp]); 

            Vector_2 s1 = query - p1;
            Vector_2 s2 = p2 - p1;
            alpha_cot.push_back(internal::cotangent_2(s2, s1, m_traits));

            s1 = p1 - p2;
            s2 = query - p2;
            beta_cot.push_back(internal::cotangent_2(s2, s1, m_traits));
          }

          FT W = FT(0);
          for (std::size_t j = 0; j < nn; ++j) {
            const std::size_t jp  = (j + 1) % nn;
            const std::size_t idx = neighbors[jp];

            const FT w = -( alpha_cot[j] + beta_cot[jp] );
            W -= w;

            if (m_domain.is_on_boundary(idx)) {
              for (std::size_t k = 0; k < n; ++k)
                b(indices[i], k) -= boundary(indices[idx], k) * w;
            } else {
              triplet_list.push_back(
                TripletFT(indices[i], indices[idx], w));
            }
          }
          triplet_list.push_back(
            TripletFT(indices[i], indices[i], W));
        }
      }

      A.setFromTriplets(
        triplet_list.begin(), triplet_list.end());
      A.makeCompressed();
      solve_linear_system(A, b, x);

      for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t i = 0; i < N; ++i) {
          if (m_domain.is_on_boundary(i))
            m_coordinates[i][k] = boundary(indices[i], k);
          else             
            m_coordinates[i][k] = x(indices[i], k);
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
      m_coordinates.clear(); m_b.clear(); m_element.clear();
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
    const Domain& m_domain;

    Solver& m_solver;

    const VertexMap m_vertex_map;
    const GeomTraits m_traits;

    std::vector<Point_2> m_polygon;
    std::vector< std::vector<FT> > m_coordinates;

    std::vector<FT> m_b;
    std::vector<std::size_t> m_element;

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
