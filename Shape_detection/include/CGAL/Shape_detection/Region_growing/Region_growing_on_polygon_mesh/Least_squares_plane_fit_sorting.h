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
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utilities.h>
#include <CGAL/Shape_detection/Region_growing/internal/property_maps.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Best least squares plane fit sorting of faces in a polygon mesh.

    This class allows to sort indices of faces in a polygon mesh, where 
    the sorting is based on the quality of the local best least squares plane fit. 
    The plane is fitted to the face vertices and vertices of its neighbors 
    found via the `Connectivity` class. The faces are sorted in the decreasing 
    quality order, that is the first index - the best quality.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam FaceListGraph 
    is a model of `FaceListGraph`.

    \tparam Connectivity 
    is a model of `RegionGrowingConnectivity`.

    \tparam FaceRange 
    is a model of `ConstRange`, whose iterator type is `RandomAccessIterator` 
    and value type is a face type used in `FaceListGraph`.

    \tparam VertexToPointMap 
    is an `LvaluePropertyMap` that maps a polygon mesh vertex to `CGAL::Point_3`.
  */
  template<
  typename GeomTraits, 
  typename FaceListGraph,
  typename Connectivity,
  typename FaceRange = typename FaceListGraph::Face_range,
  typename VertexToPointMap = typename boost::property_map<FaceListGraph, CGAL::vertex_point_t>::type>
  class Polygon_mesh_least_squares_plane_fit_sorting {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Face_graph = FaceListGraph;
    using Input_connectivity = Connectivity;
    using Face_range = FaceRange;
    using Vertex_to_point_map = VertexToPointMap;
    /// \endcond

    /// Property map that returns the quality ordered seed indices of the faces.
    using Seed_map = internal::Seed_property_map;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param polygon_mesh 
      An instance of a `FaceListGraph` that represents a polygon mesh.

      \param connectivity 
      An instance of the `Connectivity` class that is used
      internally to access face neighbors.

      \param vertex_to_point_map 
      An instance of an `LvaluePropertyMap` that maps a polygon mesh 
      vertex to `CGAL::Point_3`.

      \pre `total_number_of_faces > 0`
    */
    Polygon_mesh_least_squares_plane_fit_sorting(
      const FaceListGraph& polygon_mesh,
      Connectivity& connectivity,
      const VertexToPointMap vertex_to_point_map = VertexToPointMap()) :
    m_face_graph(polygon_mesh),
    m_connectivity(connectivity),
    m_face_range(CGAL::faces(m_face_graph)),
    m_vertex_to_point_map(vertex_to_point_map) {

      CGAL_precondition(m_face_range.size() > 0);

      m_order.resize(m_face_range.size());
      for (std::size_t i = 0; i < m_face_range.size(); ++i)
        m_order[i] = i;
      m_scores.resize(m_face_range.size());
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief Sorts face indices.
    */
    void sort() {
      
      compute_scores();
      CGAL_postcondition(m_scores.size() > 0);

      Compare_scores cmp(m_scores);
      std::sort(m_order.begin(), m_order.end(), cmp);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief Returns the `Seed_map` that allows to access 
      the ordered face indices.
    */
    Seed_map seed_map() {
      return Seed_map(m_order);
    }

    /// @}

  private:

    // Types.
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;
    using Compare_scores = internal::Compare_scores<Local_FT>;

    // Functions.
    void compute_scores() {

      std::vector<std::size_t> neighbors;
      for (std::size_t i = 0; i < m_face_range.size(); ++i) {
        
        neighbors.clear();
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
                            
            const auto& tmp_point = get(m_vertex_to_point_map, *vertex);
            points.push_back(m_to_local_converter(tmp_point));
          }
        }
        CGAL_postcondition(points.size() > 0);

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
    Input_connectivity& m_connectivity;
    const Face_range m_face_range;
    const Vertex_to_point_map m_vertex_to_point_map;

    std::vector<std::size_t> m_order;
    std::vector<Local_FT> m_scores;

    const To_local_converter m_to_local_converter;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H
