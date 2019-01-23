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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ADJACENT_FACES_CONNECTIVITY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ADJACENT_FACES_CONNECTIVITY_H

// #include <CGAL/license/Shape_detection.h>

// Boost includes.
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Face graph includes.
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_maps.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh
    \brief Find all faces that share an edge with a given face.
    \tparam FaceGraph General face graph. Model of `FaceListGraph`.
    \tparam FaceRange A random access range with graph faces.
    \cgalModels `RegionGrowingConnectivity`
  */
  template<
  typename FaceGraph, 
  typename FaceRange = typename FaceGraph::Face_range>
  class Polygon_mesh_adjacent_faces_connectivity {

  public:

    /// \name Types
    /// @{

    using Face_graph = FaceGraph;
    ///< General face graph. Model of `FaceGraph`.

    using Face_range = FaceRange;
    ///< A random access range with graph faces.

    ///< \cond SKIP_IN_MANUAL
    using Face_to_index_map 
    = internal::Item_to_index_property_map<Face_range>;
    ///< \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      The constructor that takes an arbitrary face graph.
    */
    Polygon_mesh_adjacent_faces_connectivity(
      const Face_graph& face_graph) :
    m_face_graph(face_graph),
    m_face_range(CGAL::faces(m_face_graph)),
    m_face_to_index_map(m_face_range)
    { }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      Using a query index `query_index`, this function retrieves indices 
      of all neighboring faces and stores them in `neighbors`.
      \tparam OutputIterator
    */
    template<typename OutputIterator>
    void get_neighbors(
      const std::size_t query_index, 
      OutputIterator neighbors) const {
                    
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_face_range.size());

      const auto& query_face = *(m_face_range.begin() + query_index);
      const auto& query_halfedge = CGAL::halfedge(query_face, m_face_graph);

      const auto& faces = CGAL::faces_around_face(query_halfedge, m_face_graph);
      
      for (auto face = faces.begin(); face != faces.end(); ++face) {
        const std::size_t face_index = get(m_face_to_index_map, *face);
        
        if (face_index != std::size_t(-1)) // not a null face
          *(neighbors++) = face_index;
      }
    }

    /// @}

  private:

    // Fields.
    const Face_graph& m_face_graph;
    const Face_range m_face_range;

    const Face_to_index_map m_face_to_index_map;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ADJACENT_FACES_CONNECTIVITY_H
