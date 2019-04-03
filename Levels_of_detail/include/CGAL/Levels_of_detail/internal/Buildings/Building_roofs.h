// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <memory>
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spacial search.
#include <CGAL/Levels_of_detail/internal/Spacial_search/Sphere_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_plane_fit_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_plane_fit_sorting.h>

// Partitioning.
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_3.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_3.h>

// Graphcut.
#include <CGAL/Levels_of_detail/internal/Graphcut/Graphcut.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_ground_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_roofs_estimator.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Building_roofs {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    
    using Points_3 = std::vector<std::size_t>;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using Indexer = internal::Indexer<Point_3>;

    using Vectors_3 = std::vector<Vector_3>;
    using Pair_item_3 = std::pair<Point_3, Vector_3>;
    using Pair_range_3 = std::vector<Pair_item_3>;
    using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_3>;
    using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_3>;

    using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Points_3, Point_map_3>;
    using Normal_estimator_3 = 
    internal::Estimate_normals_3<Traits, Points_3, Point_map_3, Sphere_neighbor_query>;
    using LSPF_region = 
    internal::Least_squares_plane_fit_region<Traits, Pair_range_3, First_of_pair_map, Second_of_pair_map>;
    using LSPF_sorting =
    internal::Least_squares_plane_fit_sorting<Traits, Points_3, Sphere_neighbor_query, Point_map_3>;
    using Region_growing_3 = 
    internal::Region_growing<Points_3, Sphere_neighbor_query, LSPF_region, typename LSPF_sorting::Seed_map>;

    using Building = internal::Building<Traits>;
    using Triangulation = typename Building::Base::Triangulation::Delaunay;
    using Building_ground_estimator = internal::Building_ground_estimator<Traits, Triangulation>;
    using Approximate_face = internal::Partition_edge_3<Traits>;
    using Building_walls_estimator = internal::Building_walls_estimator<Traits>;
    using Building_roofs_estimator = internal::Building_roofs_estimator<Traits, Points_3, Point_map_3>;

    /*
    using Partition_3 = internal::Partition_3<Traits>;
    using Kinetic_partitioning_3 = internal::Kinetic_partitioning_3<Traits>;
    using Visibility_3 = internal::Visibility_3<Traits>;
    using Graphcut_2 = internal::Graphcut<Traits, Partition_3>;
    using Partition_faces_3 = std::vector<typename Partition_3::Face>; */
    
    Building_roofs(
      const Data_structure& data,
      const Building& building,
      const Points_3& cluster) : 
    m_data(data),
    m_building(building),
    m_cluster(cluster),
    m_empty(false) { 
      if (cluster.empty())
        m_empty = true;
    }

    void detect_roofs() {
      if (empty())
        return;

      extract_roof_regions_3(
        m_data.parameters.buildings.region_growing_scale_3,
        m_data.parameters.buildings.region_growing_noise_level_3,
        m_data.parameters.buildings.region_growing_angle_3,
        m_data.parameters.buildings.region_growing_min_area_3,
        m_data.parameters.buildings.region_growing_distance_to_line_3,
        m_data.parameters.buildings.alpha_shape_size_2);
      make_approximate_bounds();
    }

    void compute_roofs() {
      if (empty())
        return;

    }

    const bool empty() const {
      return m_empty;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> 
    get_roof_points(
      OutputIterator output,
      std::size_t& roof_index) const {

      if (m_roof_points_3.empty())
        return boost::none;

      for (std::size_t i = 0; i < m_roof_points_3.size(); ++i) {
        for (const std::size_t idx : m_roof_points_3[i]) {
          const Point_3& p = get(m_data.point_map_3, *(m_cluster.begin() + idx));
          *(output++) = std::make_pair(p, roof_index);
        }
        ++roof_index;
      }
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_approximate_bounds(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& building_index) const {
      
      m_building_ground.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      for (const auto& wall : m_building_walls)
        wall.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      for (const auto& roof : m_building_roofs)
        roof.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

  private:
    const Data_structure& m_data;
    const Building& m_building;
    const Points_3& m_cluster;
    
    bool m_empty;
    std::vector< std::vector<std::size_t> > m_roof_points_3;
    Approximate_face m_building_ground;
    std::vector<Approximate_face> m_building_walls;
    std::vector<Approximate_face> m_building_roofs;

    void extract_roof_regions_3(
      const FT region_growing_scale_3,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3,
      const FT alpha_shape_size_2) {
        
      if (empty()) return;
      m_roof_points_3.clear();

      Sphere_neighbor_query neighbor_query(
        m_cluster, region_growing_scale_3, m_data.point_map_3);

      Vectors_3 normals;
      Normal_estimator_3 estimator(
        m_cluster, neighbor_query, m_data.point_map_3);
      estimator.get_normals(normals);

      CGAL_assertion(m_cluster.size() == normals.size());
      Pair_range_3 range;
      range.reserve(m_cluster.size());
      for (std::size_t i = 0; i < m_cluster.size(); ++i) {
        const Point_3& p = get(m_data.point_map_3, m_cluster[i]);
        range.push_back(std::make_pair(p, normals[i]));
      }

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSPF_region region(
        range, 
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        point_map,
        normal_map);

      LSPF_sorting sorting(
        m_cluster, neighbor_query, m_data.point_map_3);
      sorting.sort();

      Region_growing_3 region_growing(
        m_cluster,
        neighbor_query,
        region,
        sorting.seed_map());
      region_growing.detect(std::back_inserter(m_roof_points_3));
    }

    void make_approximate_bounds() {
        
      const FT bottom_z = m_building.bottom_z;
      const Building_ground_estimator gestimator(
        m_building.base1.triangulation.delaunay,
        bottom_z);
      gestimator.estimate(m_building_ground);

      FT top_z = m_building.top_z;
      CGAL_assertion(top_z > bottom_z);
      top_z -= (top_z - bottom_z) / FT(2);

      const Building_walls_estimator westimator(
        m_building.edges1,
        bottom_z, 
        top_z);
      westimator.estimate(m_building_walls);

      const Building_roofs_estimator restimator(
        m_cluster,
        m_data.point_map_3,
        m_roof_points_3);
      restimator.estimate(m_building_roofs);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H
