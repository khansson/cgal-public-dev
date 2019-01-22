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
    \ingroup PkgShapeDetectionRGOnGraph
    \brief Local and global conditions for the region growing algorithm on a face graph.
    \tparam Traits Model of `Kernel`
    \tparam FaceGraph General face graph. Model of `FaceGraph`.
    \tparam FaceRange An arbitrary range with graph faces, given an IndexToFaceMap is provided.
    \tparam VertexToPointMap An `LvaluePropertyMap` that maps a graph vertex to `Point_3`.
    \cgalModels `RegionGrowingPropagationConditions`
  */
  template<
  typename GeomTraits, 
  typename FaceGraph,
  typename FaceRange = typename FaceGraph::Face_range,
  typename VertexToPointMap = typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type>
  class Polygon_mesh_least_squares_plane_fit_conditions {

  public:

    /// \name Types
    /// @{

    using Traits = GeomTraits;

    using Face_graph = FaceGraph;
    ///< General face graph. Model of `FaceGraph`.

    using Face_range = FaceRange;
    ///< An arbitrary range with graph faces.

    using Vertex_to_point_map = VertexToPointMap;
    ///< An `LvaluePropertyMap` that maps a graph vertex to `Point_3`.

    using FT = typename Traits::FT; ///< Number type
    using Point_3 = typename Traits::Point_3; ///< Point type
    using Vector_3 = typename Traits::Vector_3; ///< Vector type
    using Plane_3 = typename Traits::Plane_3; ///< Plane type

    ///< \cond SKIP_IN_MANUAL
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
    ///< \endcond

    /// @}

    /*!
      Each region is represented by a plane. 
      The constructor requires an input range with graph faces 
      and three parameters can be provided, in order: 
      the maximum distance from a point to the region, 
      the minimum dot product between the normal associated with the point and the normal of the region, 
      and the minimum number of points a region must have.
      In addition, you can provide instances of the Vertex_to_point_map and Traits classes.
    */
    Polygon_mesh_least_squares_plane_fit_conditions(
      const Face_graph& face_graph,
      const FT distance_threshold = FT(1), 
      const FT normal_threshold = FT(9) / FT(10), 
      const std::size_t min_region_size = 1, 
      const Vertex_to_point_map vertex_to_point_map = Vertex_to_point_map(), 
      const Traits traits = Traits()) :
    m_face_graph(face_graph),
    m_face_range(CGAL::faces(m_face_graph)),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(normal_threshold),
    m_min_region_size(min_region_size),
    m_vertex_to_point_map(vertex_to_point_map),
    m_squared_length_3(traits.compute_squared_length_3_object()),
    m_squared_distance_3(traits.compute_squared_distance_3_object()),
    m_scalar_product_3(traits.compute_scalar_product_3_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(distance_threshold >= FT(0));
      CGAL_precondition(normal_threshold >= FT(0) && normal_threshold <= FT(1));
      CGAL_precondition(min_region_size > 0);
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      Local conditions that check if a query item belongs to the given region.
      \tparam Region CGAL::Shape_detection::Region_growing::Region
    */
    template<typename ItemRange>
    bool belongs_to_region(
      const std::size_t query_index, 
      const ItemRange& region) const {

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
      Global conditions that check if a region size is large enough to be accepted.
      \tparam Region CGAL::Shape_detection::Region_growing::Region
    */
    template<typename ItemRange>
    inline bool is_valid_region(const ItemRange& region) const {
      return ( region.size() >= m_min_region_size );
    }

    /*!
      Update the class's best fit plane that will be used later by local conditions.
      \tparam Region CGAL::Shape_detection::Region_growing::Region
    */
    template<typename ItemRange>
    void update(const ItemRange& region) {

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

        std::vector<Local_point_3> points(region.size());
        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_face_range.size());

          const auto& face = *(m_face_range.begin() + region[i]);
          const auto& halfedge = CGAL::halfedge(face, m_face_graph);

          const auto& vertices = CGAL::vertices_around_face(halfedge, m_face_graph);
          for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex) {
                            
            const Point_3& tmp_point = get(m_vertex_to_point_map, *vertex);
            points[i] = m_to_local_converter(tmp_point);
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
