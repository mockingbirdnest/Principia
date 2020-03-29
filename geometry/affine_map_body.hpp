
#pragma once

#include <utility>

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"
#include "affine_map.hpp"

namespace principia {
namespace geometry {
namespace internal_affine_map {

// The map is represented as x ↦ linear_map(x - from_origin) + to_origin.  This
// numerically better behaved than x ↦ linear_map(x) + translation with
// translation = to_origin - linear_map(from_origin).
template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::AffineMap(
    Point<FromVector> const& from_origin,
    Point<ToVector> const& to_origin,
    LinearMap<FromFrame, ToFrame> linear_map)
    : from_origin_(from_origin),
      to_origin_(to_origin),
      linear_map_(std::move(linear_map)) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<ToFrame, FromFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::Inverse() const {
  return AffineMap<ToFrame, FromFrame, Scalar, LinearMap>(
      /*from_origin=*/to_origin_,
      /*to_origin=*/from_origin_,
      /*linear_map=*/linear_map_.Inverse());
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

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
template<typename F, typename T, typename>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::Identity() {
  return AffineMap(Point<FromVector>(),
                   Point<ToVector>(),
                   LinearMap<FromFrame, ToFrame>::Identity());
}

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
LinearMap<FromFrame, ToFrame> const&
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::linear_map() const {
  return linear_map_;
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
void AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::WriteToMessage(
    not_null<serialization::AffineMap*> const message) const {
  FromFrame::WriteToMessage(message->mutable_from_frame());
  ToFrame::WriteToMessage(message->mutable_to_frame());
  from_origin_.WriteToMessage(message->mutable_from_origin());
  to_origin_.WriteToMessage(message->mutable_to_origin());
  linear_map_.WriteToMessage(message->mutable_linear_map());
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
template<typename, typename, typename>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ReadFromMessage(
    serialization::AffineMap const& message) {
  FromFrame::ReadFromMessage(message.from_frame());
  ToFrame::ReadFromMessage(message.to_frame());
  return AffineMap(Point<FromVector>::ReadFromMessage(message.from_origin()),
                   Point<ToVector>::ReadFromMessage(message.to_origin()),
                   LinearMap<FromFrame, ToFrame>::ReadFromMessage(
                       message.linear_map()));
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ThroughFrame, Scalar, LinearMap> const& right) {
  return AffineMap<FromFrame, ToFrame, Scalar, LinearMap>(
      /*from_origin=*/right.from_origin_,
      /*to_origin=*/left.to_origin_ +
          left.linear_map_(right.to_origin_ - left.from_origin_),
      /*linear_map=*/left.linear_map_ * right.linear_map_);
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
std::ostream& operator<<(
    std::ostream& out,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& affine_map) {
  return out << "{from: " << affine_map.from_origin_
             << ", to: " << affine_map.to_origin_
             << ", map: " << affine_map.linear_map_ << "}";
}

}  // namespace internal_affine_map
}  // namespace geometry
}  // namespace principia
