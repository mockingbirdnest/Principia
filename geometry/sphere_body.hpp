
#pragma once

#include "geometry/sphere.hpp"

namespace principia {
namespace geometry {
namespace internal_sphere {

template<typename Scalar, typename Frame>
Sphere<Scalar, Frame>::Sphere(Point<Vector<Scalar, Frame>> const& centre,
                              Scalar const& radius)
    : centre_(centre),
      radius_(radius),
      radius²_(radius_ * radius_) {}

template<typename Scalar, typename Frame>
Point<Vector<Scalar, Frame>> const& Sphere<Scalar, Frame>::centre() const {
  return centre_;
}

template<typename Scalar, typename Frame>
Scalar const& Sphere<Scalar, Frame>::radius() const {
  return radius_;
}

template<typename Scalar, typename Frame>
Product<Scalar, Scalar> const Sphere<Scalar, Frame>::radius²() const {
  return radius²_;
}

}  // namespace internal_sphere
}  // namespace geometry
}  // namespace principia
