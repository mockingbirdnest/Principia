
#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/rp2_element.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

// A perspective using the pinhole camera model.  It project a point of
// |FromFrame| to an element of ℝP².  |ToFrame| is the frame of the camera.  In
// that frame the camera is located at the origin and looking at the positive
// z-axis.  The x- and y- axis of the camera correspond to those of ℝP².
template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
class Perspective final {
 public:
  Perspective(
      AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& to_camera,
      Scalar const& focal);
  Perspective(
      AffineMap<ToFrame, FromFrame, Scalar, LinearMap> const& from_camera,
      Scalar const& focal);

  RP2Element<Scalar> operator()(
      Point<Vector<Scalar, FromFrame>> const& point) const;

 private:
  AffineMap<FromFrame, ToFrame, Scalar, LinearMap> to_camera_;
  Scalar focal_;
};

}  // namespace internal_perspective

using internal_perspective::Perspective;

}  // namespace geometry
}  // namespace principia

#include "geometry/perspective_body.hpp"
