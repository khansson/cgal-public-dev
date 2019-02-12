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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_CONDITIONS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_CONDITIONS_H

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

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Least squares plane fit conditions on a polygon mesh.

    This class implements propagation conditions for detecting planes 
    on a polygon mesh via the `Shape_detection::Region_growing` approach, 
    where quality of detected planes is based on a local least squares plane fit.

    \tparam GeomTraits 
    is a model of `Kernel`.

    \tparam FaceListGraph 
    is a model of `FaceListGraph`.

    \tparam FaceRange 
    is a model of `ConstRange`, whose iterator type is `RandomAccessIterator` 
    and value type is a face type used in `FaceListGraph`.

    \tparam VertexToPointMap 
    is an `LvaluePropertyMap` that maps a polygon mesh vertex to `CGAL::Point_3`.
    
    \cgalModels `RegionGrowingPropagationConditions`
  */
  template<
  typename GeomTraits, 
  typename FaceListGraph,
  typename FaceRange = typename FaceListGraph::Face_range,
  typename VertexToPointMap = typename boost::property_map<FaceListGraph, CGAL::vertex_point_t>::type>
  class Polygon_mesh_least_squares_plane_fit_conditions {

  public:

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Face_graph = FaceListGraph;
    using Face_range = FaceRange;
    using Vertex_to_point_map = VertexToPointMap;

    using Vector_3 = typename Traits::Vector_3;

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_3 = typename Traits::Compute_squared_length_3;
    using Squared_distance_3 = typename Traits::Compute_squared_distance_3;
    using Scalar_product_3 = typename Traits::Compute_scalar_product_3;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;
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

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes all internal data structures.

      \param polygon_mesh 
      An instance of a `FaceListGraph` that represents
      a polygon mesh.

      \param distance_threshold 
      Maximum distance from a face to the region represented by a plane of 
      the type `CGAL::Plane_3`.

      \param normal_threshold 
      Minimum dot product between the face normal and the normal assigned 
      to the region represented by a plane of the type `CGAL::Plane_3`.

      \param min_region_size 
      The minimum number of faces a region must have.
      
      \param vertex_to_point_map 
      An instance of an `LvaluePropertyMap` that maps a polygon mesh 
      vertex to `CGAL::Point_3`.

      \param traits
      An instance of the `GeomTraits` class.

      \pre `distance_threshold >= 0`
      \pre `normal_threshold >= 0 && normal_threshold <= 1`
      \pre `min_region_size > 0`
    */
    Polygon_mesh_least_squares_plane_fit_conditions(
      const FaceListGraph& polygon_mesh,
      const FT distance_threshold = FT(1), 
      const FT normal_threshold = FT(9) / FT(10), 
      const std::size_t min_region_size = 1, 
      const VertexToPointMap vertex_to_point_map = VertexToPointMap(), 
      const GeomTraits traits = GeomTraits()) :
    m_face_graph(polygon_mesh),
    m_face_range(CGAL::faces(m_face_graph)),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(normal_threshold),
    m_min_region_size(min_region_size),
    m_vertex_to_point_map(vertex_to_point_map),
    m_squared_length_3(traits.compute_squared_length_3_object()),
    m_squared_distance_3(traits.compute_squared_distance_3_object()),
    m_scalar_product_3(traits.compute_scalar_product_3_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(m_face_range.size() > 0);

      CGAL_precondition(m_distance_threshold >= FT(0));
      CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));
      CGAL_precondition(m_min_region_size > 0);
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief Checks if a face belongs to a region.

      Checks if the face with the index `query_index` belongs to the region
      that is currently getting developed using the `distance_threshold` and
      `normal_threshold` values.

      \param query_index
      Index of the query face.

      The second parameter is not used in this implementation.

      Implements the function `RegionGrowingPropagationConditions::belongs_to_region()`.

      \pre `query_index >= 0 && query_index < total_number_of_faces`
    */
    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_face_range.size());

      const auto& face = *(m_face_range.begin() + query_index);
                
      Vector_3 face_normal;
      get_face_normal(face, face_normal);

      const FT distance_to_fitted_plane = 
      get_max_face_distance(face);
      
      const FT cos_angle = 
      CGAL::abs(m_scalar_product_3(face_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_plane <= m_distance_threshold ) && 
        ( cos_angle >= m_normal_threshold ));
    }

    /*!
      \brief Validates the final `region`.

      Controls if the `region` that has been created contains at least 
      `min_region_size` faces.

      \param region
      Stores indices of all faces that belong to the region.

      Implements the function `RegionGrowingPropagationConditions::is_valid_region()`.
    */
    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return ( region.size() >= m_min_region_size );
    }

    /*!
      \brief Recomputes a least squares plane.

      Recomputes the internal least squares plane that represents the `region`
      currently being developed. The plane is fitted to the polygon mesh
      vertices of all faces, which have been added to the `region` so far.

      \param region
      Stores indices of all faces that belong to the region.

      Implements the function `RegionGrowingPropagationConditions::update()`.

      \pre `region.size() > 0`
    */
    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference plane and normal

        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_face_range.size());

        const auto& face = *(m_face_range.begin() + region[0]);
        Point_3 face_centroid;

        get_face_centroid(face, face_centroid);
        get_face_normal(face, m_normal_of_best_fit);

        m_plane_of_best_fit = 
        Plane_3(face_centroid, m_normal_of_best_fit);

      } else {

        std::vector<Local_point_3> points;
        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_face_range.size());

          const auto& face = *(m_face_range.begin() + region[i]);
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
          linear_least_squares_fitting_3(
            points.begin(), points.end(), 
            fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 3>());
        #else 
          linear_least_squares_fitting_3(
            points.begin(), points.end(), 
            fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 3>());
        #endif

        m_plane_of_best_fit = 
        Plane_3(
          static_cast<FT>(fitted_plane.a()), 
          static_cast<FT>(fitted_plane.b()), 
          static_cast<FT>(fitted_plane.c()), 
          static_cast<FT>(fitted_plane.d()));

        const Vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
        const FT normal_length = m_sqrt(m_squared_length_3(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

    /// @}

  private:

    template<typename Face>  
    void get_face_centroid(
      const Face& face, 
      Point_3& face_centroid) const {

      const auto& halfedge = CGAL::halfedge(face, m_face_graph);
      const auto& vertices = CGAL::vertices_around_face(halfedge, m_face_graph);

      // Compute centroid.
      FT sum = FT(0);

      FT x = FT(0);
      FT y = FT(0);
      FT z = FT(0);
                
      for (auto vertex = vertices.begin(); 
      vertex != vertices.end(); 
      ++vertex, sum += FT(1)) {
        
        const Point_3& point = get(m_vertex_to_point_map, *vertex);

        x += point.x();
        y += point.y();
        z += point.z();
      }
      CGAL_precondition(sum > FT(0));

      x /= sum;
      y /= sum;
      z /= sum;

      face_centroid = Point_3(x, y, z);
    }

    template<typename Face>
    void get_face_normal(
      const Face& face, 
      Vector_3& face_normal) const {

      const auto& halfedge = CGAL::halfedge(face, m_face_graph);
      const auto& vertices = CGAL::vertices_around_face(halfedge, m_face_graph);

      auto vertex = vertices.begin();
      const Point_3& point1 = get(m_vertex_to_point_map, *vertex); ++vertex;
      const Point_3& point2 = get(m_vertex_to_point_map, *vertex); ++vertex;
      const Point_3& point3 = get(m_vertex_to_point_map, *vertex);

      const Vector_3 tmp_normal = CGAL::normal(point1, point2, point3);
      const FT tmp_normal_length = m_sqrt(m_squared_length_3(tmp_normal));

      CGAL_precondition(tmp_normal_length > FT(0));
      face_normal = tmp_normal / tmp_normal_length;
    }

    template<typename Face>
    FT get_max_face_distance(const Face& face) const {

      const auto& halfedge = CGAL::halfedge(face, m_face_graph);
      const auto& vertices = CGAL::vertices_around_face(halfedge, m_face_graph);

      FT max_face_distance = -FT(1);
      for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex) {
        
        const Point_3& point = 
        get(m_vertex_to_point_map, *vertex);

        const FT distance = 
        m_sqrt(m_squared_distance_3(point, m_plane_of_best_fit));
        
        max_face_distance = 
        CGAL::max(distance, max_face_distance);
      }
      CGAL_postcondition(max_face_distance != -FT(1));

      return max_face_distance;
    }

    // Fields.
    const Face_graph& m_face_graph;
    const Face_range m_face_range;

    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const std::size_t m_min_region_size;
            
    const Vertex_to_point_map m_vertex_to_point_map;

    const Squared_length_3 m_squared_length_3;
    const Squared_distance_3 m_squared_distance_3;
    const Scalar_product_3 m_scalar_product_3;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Plane_3 m_plane_of_best_fit;
    Vector_3 m_normal_of_best_fit;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_CONDITIONS_H
