#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H

// #include <CGAL/license/Shape_regularization.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<typename Point_2>
  Point_2 compute_middle_point(const Point_2& source, const Point_2& target) {
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    const FT half = FT(1) / FT(2);
    const FT x = half * (source.x() + target.x());
    const FT y = half * (source.y() + target.y());
    return Point_2(x, y);
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::Vector_2
  compute_direction(const Segment_2& segment) {
    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector = typename Traits::Vector_2;

    Vector v = segment.to_vector(); 
    if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0))) 
      v = -v;
    return v;
  }
  
  template<typename Vector>
  typename Kernel_traits<Vector>::Kernel::FT
  compute_orientation(const Vector& v) {
    using Traits = typename Kernel_traits<Vector>::Kernel;
    using FT = typename Traits::FT;

    const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
    FT orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
    if (orientation < FT(0)) 
      orientation += FT(180);
    return orientation;
  }

  template<typename Point_2, typename FT>
  Point_2 transform_coordinates(const Point_2 & barycentre, const Point_2 & frame_origin, const FT angle) {

    const FT cos_val = static_cast<FT>(cos(CGAL_PI * CGAL::to_double(angle) / 180.0));
    const FT sin_val = static_cast<FT>(sin(CGAL_PI * CGAL::to_double(angle) / 180.0));

    const FT x = (barycentre.x() - frame_origin.x()) * cos_val + (barycentre.y() - frame_origin.y()) * sin_val;
    const FT y = (barycentre.y() - frame_origin.y()) * cos_val - (barycentre.x() - frame_origin.x()) * sin_val;

    return Point_2(x, y);
  }


} // internal
} // Regularization
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
