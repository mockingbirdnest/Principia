
#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/rp2_element.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

//TODO(phl):comments
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
