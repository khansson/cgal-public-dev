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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H

#include <CGAL/license/Shape_detection.h>

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
namespace Polygon_mesh {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Adjacent faces connectivity on a polygon mesh.

    This class finds all faces, which are adjacent to a given face on a 
    polygon mesh and can be used to parameterize the `Shape_detection::Region_growing` 
    algorithm. The polygon mesh is represented as a `FaceListGraph`.

    \tparam FaceListGraph 
    is a model of `FaceListGraph`.

    \tparam FaceRange 
    is a model of `ConstRange`, whose iterator type is `RandomAccessIterator` 
    and value type is a face type used in `FaceListGraph`.

    \cgalModels `RegionGrowingConnectivity`
  */
  template<
  typename FaceListGraph, 
  typename FaceRange = typename FaceListGraph::Face_range>
  class One_ring_neighbor_query {

  public:

    /// \cond SKIP_IN_MANUAL
    using Face_graph = FaceListGraph;
    using Face_range = FaceRange;

    using Face_to_index_map 
    = internal::Item_to_index_property_map<Face_range>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param polygon_mesh An instance of a `FaceListGraph` that represents
      a polygon mesh.

      \pre `total_number_of_faces > 0`
    */
    One_ring_neighbor_query(
      const FaceListGraph& polygon_mesh) :
    m_face_graph(polygon_mesh),
    m_face_range(CGAL::faces(m_face_graph)),
    m_face_to_index_map(m_face_range) { 

      CGAL_precondition(m_face_range.size() > 0);
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief Returns adjacent faces to a given face.

      This function returns indices of all faces, 
      which are adjacent to the face with the index `query_index`. 
      These indices are returned in `neighbors`.

      \param query_index
      Index of the query face.

      \param neighbors
      An `std::vector<std::size_t>` with the indices of faces, which are 
      adjacent to the face with the index `query_index`.

      Implements the function `RegionGrowingConnectivity::neighbors()`.

      \pre `query_index >= 0 && query_index < total_number_of_faces`
    */
    void operator()(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
                    
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_face_range.size());

      const auto& query_face = *(m_face_range.begin() + query_index);
      const auto& query_halfedge = CGAL::halfedge(query_face, m_face_graph);

      const auto& faces = CGAL::faces_around_face(query_halfedge, m_face_graph);
      for (const auto& face : faces) {
        const std::size_t face_index = get(m_face_to_index_map, face);
        
        if (face_index != std::size_t(-1)) // not a null face
          neighbors.push_back(face_index);
      }
    }

    /// @}

  private:

    // Fields.
    const Face_graph& m_face_graph;
    const Face_range m_face_range;

    const Face_to_index_map m_face_to_index_map;
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H
