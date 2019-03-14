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

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_weights_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  template<
  typename VertexRange,
  typename GeomTraits,
  typename InputDomain,
  typename InputSolver,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Harmonic_coordinates_2 : public Discrete_harmonic_weights_2<VertexRange, GeomTraits, VertexMap> {

    public:
      using Base = Discrete_harmonic_weights_2<VertexRange, GeomTraits, VertexMap>;

      using Vertex_range = VertexRange;
      using Traits = GeomTraits;
      using Domain = InputDomain;
      using Solver = InputSolver;
      using Vertex_map = VertexMap;
      
      using Point_2 = typename GeomTraits::Point_2;

      Harmonic_coordinates_2() {

      }

      template<typename OutputIterator>
      boost::optional<OutputIterator> operator()(
        const Point_2& p, 
        OutputIterator weights) const {

        // evaluate
      }

      bool is_valid_point(const Point_2& query) const {

        // override
      }

      void compute() {

      }

      template<typename OutputIterator>
      boost::optional<OutputIterator> coordinates(OutputIterator output) {

      }

      template<typename OutputIterator>
      boost::optional<OutputIterator> coordinates(
        const std::size_t vertex_index, 
        OutputIterator output) {

      }

      void clear() {

      }

      void release_memory() {

      }

    private:
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_HARMONIC_COORDINATES_2_H
