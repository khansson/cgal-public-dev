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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H

// #include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// Boost includes.
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Face graph includes.
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

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

  template<
  typename GeomTraits, 
  typename FaceListGraph,
  typename Connectivity_,
  typename FaceRange = typename FaceListGraph::Face_range,
  typename VertexToPointMap = typename boost::property_map<FaceListGraph, CGAL::vertex_point_t>::type>
  class Polygon_mesh_least_squares_plane_fit_sorting {

  public:

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Face_graph = FaceListGraph;
    using Face_range = FaceRange;
    using Vertex_to_point_map = VertexToPointMap;
    using Input_connectivity = Connectivity_;

    /// \endcond

    /// \name Types
    /// @{

    /// Number type.
    using FT = typename GeomTraits::FT;

    /// Point type.
    using Point_3 = typename GeomTraits::Point_3; 

    /// Type of the plane.
    using Plane_3 = typename GeomTraits::Plane_3;

    /// @}

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

    Polygon_mesh_least_squares_plane_fit_sorting(
      const FaceListGraph& polygon_mesh,
      Input_connectivity& connectivity,
      const VertexToPointMap vertex_to_point_map = VertexToPointMap(), 
      const GeomTraits traits = GeomTraits()) :
    m_face_graph(polygon_mesh),
    m_connectivity(connectivity),
    m_face_range(CGAL::faces(m_face_graph)),
    m_vertex_to_point_map(vertex_to_point_map) {

      CGAL_precondition(m_face_range.size() > 0);
      std::size_t input_size = m_face_range.end() - m_face_range.begin();
      m_order.resize(input_size);
      m_scores.resize(input_size);
      for (std::size_t i = 0; i < input_size; ++i) m_order[i] = i;

    }

    void sort() {
      calculate_scores();
      Compare cmp(m_scores);
      std::sort(m_order.begin(), m_order.end(), cmp);
    }

    Seed_map seed_map() {
      return Seed_map(m_order);
    }

  private:

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

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

      for (std::size_t i = 0; i < m_face_range.size(); ++i) {
        
        std::vector<std::size_t> neighbors;
        m_connectivity.neighbors(i, neighbors);
        neighbors.push_back(i);

        std::vector<Local_point_3> points;

        for (std::size_t j = 0; j < neighbors.size(); ++j) {

          CGAL_precondition(neighbors[j] >= 0);
          CGAL_precondition(neighbors[j] < m_face_range.size());

          const auto& face = *(m_face_range.begin() + neighbors[j]);
          const auto& halfedge = CGAL::halfedge(face, m_face_graph);

          const auto& vertices = CGAL::vertices_around_face(halfedge, m_face_graph);
          for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex) {
                            
            const Point_3& tmp_point = get(m_vertex_to_point_map, *vertex);
            points.push_back(m_to_local_converter(tmp_point));

          }
        }

        CGAL_precondition(points.size() > 0);

        Local_plane_3 fitted_plane;
        Local_point_3 fitted_centroid;

        #ifndef CGAL_EIGEN3_ENABLED
          m_scores[i] = linear_least_squares_fitting_3(
            points.begin(), points.end(), 
            fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 3>());
        #else 
          m_scores[i] = linear_least_squares_fitting_3(
            points.begin(), points.end(), 
            fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 3>());
        #endif
        
      }

    }

    // Fields.
    const Face_graph& m_face_graph;
    const Face_range m_face_range;
            
    const Vertex_to_point_map m_vertex_to_point_map;

    const To_local_converter m_to_local_converter;

    Input_connectivity& m_connectivity;
    std::vector<std::size_t> m_order;
    std::vector<Local_FT> m_scores;
    
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H
