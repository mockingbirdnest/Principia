
#pragma once

#include "geometry/perspective.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::Perspective(
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& to_camera,
    Scalar const& focal)
    : to_camera_(to_camera),
      focal_(focal) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::Perspective(
    AffineMap<ToFrame, FromFrame, Scalar, LinearMap> const& from_camera,
    Scalar const& focal)
    : to_camera_(from_camera.Inverse()),
      focal_(focal) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
RP2Point<Scalar, ToFrame> Perspective<FromFrame, ToFrame, Scalar, LinearMap>::
operator()(Point<Vector<Scalar, FromFrame>> const& point) const {
  Point<Vector<Scalar, ToFrame>> const point_in_camera = to_camera_(point);
  Vector<Scalar, ToFrame> const displacement_in_camera =
      point_in_camera - ToFrame::origin;
  R3Element<Scalar> const coordinates_in_camera =
      displacement_in_camera.coordinates();
  // This is the actual pinhole camera projection.
  return RP2Point<Scalar, ToFrame>(coordinates_in_camera.x,
                                   coordinates_in_camera.y,
                                   coordinates_in_camera.z / focal_);
}

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
bool Perspective<FromFrame, ToFrame, Scalar, LinearMap>::IsHiddenBySphere(
    Point<Vector<Scalar, FromFrame>> const& center,
    Scalar const& radius) const {}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
