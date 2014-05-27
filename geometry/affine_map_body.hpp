#pragma once

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {

// Since the map is represented as x -> linear_map_(x) + translation_ and since
// it is x -> linear_map(x - from_origin) + to_origin, we get
// linear_map_ = linear_map, translation_ = to_origin - linear_map(from_origin).
template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::AffineMap(
    Point<FromVector> const& from_origin,
    Point<ToVector> const& to_origin,
    LinearMap<FromFrame, ToFrame> const& linear_map)
    : linear_map_(linear_map),
      translation_(to_origin.coordinates_ -
                       linear_map(from_origin.coordinates_)) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<ToFrame, FromFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::Inverse() const {
  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> result;
  result.linear_map_ = linear_map_.Inverse();
  result.translation_ = -result.linear_map_(translation_);
  return result;
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
Point<typename AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ToVector>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::operator()(
    Point<FromVector> const& point) const {
  return Point<
      typename AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ToVector>(
          linear_map_(point.coordinates_) + translation_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& right) {
  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> result;
  result.linear_map_ = left.linear_map_ * right.linear_map_;
  result.translation_ = left.linear_map_(right.translation_) +
                        left.translation_;
  return result;
}

}  // namespace geometry
}  // namespace principia
