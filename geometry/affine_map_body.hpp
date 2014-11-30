#pragma once

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {

// The map is represented as x ↦ linear_map(x - from_origin) + to_origin.  This
// numerically better behaved than x ↦ linear_map(x) + translation with
// translation = to_origin - linear_map(from_origin).
template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::AffineMap(
    Point<FromVector> const& from_origin,
    Point<ToVector> const& to_origin,
    LinearMap<FromFrame, ToFrame> const& linear_map)
    : from_origin_(from_origin),
      to_origin_(to_origin),
      linear_map_(linear_map) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<ToFrame, FromFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::Inverse() const {
  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> result;
  result.from_origin_ = to_origin_;
  result.to_origin_ = from_origin_;
  result.linear_map_ = linear_map_.Inverse();
  return result;
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
Point<typename AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ToVector>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::operator()(
    Point<FromVector> const& point) const {
  return Point<
      typename AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ToVector>(
          linear_map_(point - from_origin_) + to_origin_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& right) {
  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> result;
  result.from_origin_ = right.from_origin_;
  result.to_origin_ = left.to_origin_ +
                      left.linear_map_(right.to_origin_ - left.from_origin_);
  result.linear_map_ = left.linear_map_ * right.linear_map_;
  return result;
}

}  // namespace geometry
}  // namespace principia
