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

#ifndef CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H
#define CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// CGAL includes.
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*! 
    \ingroup PkgBarycentric_coordinates_2DD

    \brief Delaunay domain 2.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam GeomTraits 
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap 
    is an `LvaluePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \cgalModels `DiscretizedDomain_2`
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Delaunay_domain_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL 
    using Polygon_    = Polygon;
    using GeomTraits_ = GeomTraits;
    using VertexMap_  = VertexMap;

    struct VI {
      internal::Edge_case type = internal::Edge_case::UNBOUNDED;
      std::size_t index = std::size_t(-1);
    };

    using FB  = CGAL::Delaunay_mesh_face_base_2<GeomTraits>;
    using VB  = CGAL::Triangulation_vertex_base_with_info_2<VI, GeomTraits>;
    using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<GeomTraits, TDS>;

    using Vertex_handle = typename CDT::Vertex_handle;

    using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
    using Mesher   = CGAL::Delaunay_mesher_2<CDT, Criteria>;
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

      \param polygon
      An instance of `Polygon` with vertices of a simple polygon.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon` 
      to `Point_2`.

      \param traits
      An instance of `GeomTraits`.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
    */
    Delaunay_domain_2(
      const Polygon& polygon,
      const VertexMap vertex_map = VertexMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_polygon(polygon),
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
    }

    /*!
      creates input domain.
    */
    void create(
      const FT edge_size,
      const std::list<Point_2>& list_of_seeds) {

      // Create Delaunay triangulation.
      m_cdt.clear(); m_vhs.clear();
      m_vhs.reserve(m_polygon.size());
      for (const auto& vertex : m_polygon) 
        m_vhs.push_back(m_cdt.insert(vertex));

      for(std::size_t i = 0; i < m_vhs.size(); ++i) {
        const std::size_t ip = (i + 1) % m_vhs.size();
        m_cdt.insert_constraint(m_vhs[i], m_vhs[ip]);
      }

      // Refine this triangulation.
      Mesher mesher(m_cdt);
      mesher.set_seeds(list_of_seeds.begin(), list_of_seeds.end(), true);
      mesher.set_criteria(Criteria(edge_size, edge_size));
      mesher.refine_mesh();

      // Find interior points.
      std::set<Vertex_handle> vhs;
      for (auto fh = m_cdt.finite_faces_begin(); 
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (fh->is_in_domain()) {
          vhs.insert(fh->vertex(0));
          vhs.insert(fh->vertex(1));
          vhs.insert(fh->vertex(2));
        }
      }

      m_vhs.clear(); std::size_t count = 0;
      for (auto vh : vhs) {
        const auto result = internal::locate_wrt_polygon_2(
          m_polygon, vh->point(), m_traits);
        const auto location = (*result).first;

        /*
        if (
          location == Query_point_location::ON_UNBOUNDED_SIDE)
          vh->info().type = internal::Edge_case::BOUNDARY; */

        if (
          location == Query_point_location::ON_VERTEX ||
          location == Query_point_location::ON_EDGE)
          vh->info().type = internal::Edge_case::BOUNDARY;

        if (
          location == Query_point_location::ON_BOUNDED_SIDE)
          vh->info().type = internal::Edge_case::INTERIOR;

        vh->info().index = count;
        m_vhs.push_back(vh); ++count;
      }
    }

    /*!
      get barycenters.
    */
    void get_barycenters(
      std::vector<Point_2>& barycenters) const {

      barycenters.clear();
      barycenters.reserve(m_cdt.number_of_faces());

      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        const Point_2 b = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));
        barycenters.push_back(b);
      }
    }

    /*!
      returns number of vertices.
    */
    const std::size_t number_of_vertices() const {
      return m_vhs.size();
    }

    /*!
      returns a vertex.
    */
    const Point_2& vertex(
      const std::size_t query_index) const {

      CGAL_precondition(
        query_index >= 0 && query_index < number_of_vertices());
      const auto vh = m_vhs[query_index];
      return vh->point();
    }

    /*!
      controls the boundary.
    */
    const bool is_on_boundary(
      const std::size_t query_index) const {

      CGAL_precondition(
        query_index >= 0 && query_index < number_of_vertices());
      const auto vh = m_vhs[query_index];
      return vh->info().type == internal::Edge_case::BOUNDARY;
    }

    /*!
      returns neighbors.
    */
    void operator()(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
      
      CGAL_precondition(
        query_index >= 0 && query_index < number_of_vertices());

      neighbors.clear();
      const auto vh = m_vhs[query_index];
      CGAL_assertion(vh->info().type == internal::Edge_case::INTERIOR);

      auto circ = m_cdt.incident_vertices(vh);
      const auto end = circ;
      CGAL_assertion(!circ.is_empty());

      do {
        CGAL_assertion(
          circ->info().type == internal::Edge_case::INTERIOR ||
          circ->info().type == internal::Edge_case::BOUNDARY);
        CGAL_assertion(
          !m_cdt.is_infinite(circ));

        neighbors.push_back(circ->info().index);
        ++circ;
      } while (circ != end);
      
      CGAL_assertion(neighbors.size() > 0);
    }

    /*!
      locates the point.
    */
    bool locate(
      const Point_2& query, 
      std::vector<std::size_t>& element) const {

      element.clear();
      const auto fh = m_cdt.locate(query);
      for (std::size_t i = 0; i < 3; ++i) {
        const auto vh = fh->vertex(i);
        if (
          vh->info().type == internal::Edge_case::INTERIOR ||
          vh->info().type == internal::Edge_case::BOUNDARY )
          element.push_back(vh->info().index);
      }
      return ( element.size() == 3 );
    }

  private:
      
    // Fields.
    const Polygon& m_input_polygon;
    const VertexMap m_vertex_map;
    const GeomTraits m_traits;

    std::vector<Point_2> m_polygon;
    std::vector<Vertex_handle> m_vhs;
    CDT m_cdt;
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H
