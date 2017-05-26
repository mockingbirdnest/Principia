
#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
class Perspective final {
 public:
  using FromVector = Vector<Scalar, FromFrame>;
  using ToVector = Vector<Scalar, ToFrame>;

  Perspective(
      AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& to_camera,
      Scalar const& focal);
  Perspective(
      AffineMap<ToFrame, FromFrame, Scalar, LinearMap> const& from_camera,
      Scalar const& focal);

  R3Element<Scalar> operator()(Point<FromVector> const& point) const;
};

}  // namespace geometry
}  // namespace principia
