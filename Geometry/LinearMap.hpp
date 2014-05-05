#pragma once

#include "Grassmann.hpp"

namespace principia {
namespace geometry {

template<typename Scalar,  typename FromFrame, typename ToFrame>
class LinearMap {
 public:
  LinearMap() = default;
  virtual ~LinearMap() = default;

  virtual Vector<Scalar, ToFrame> ActOn(
      Vector<Scalar, FromFrame> const& vector) const = 0;
  virtual Bivector<Scalar, ToFrame> ActOn(
      Bivector<Scalar, FromFrame> const& bivector) const = 0;
  virtual Trivector<Scalar, ToFrame> ActOn(
      Trivector<Scalar, FromFrame> const& trivector) const = 0;
};

}  // namespace geometry
}  // namespace principia
