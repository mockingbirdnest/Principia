
#pragma once

#include "geometry/sphere.hpp"

namespace principia {
namespace geometry {
namespace internal_sphere {

template<typename Frame>
Sphere<Frame>::Sphere(Position<Frame> const& centre,
                      Length const& radius)
    : centre_(centre),
      radius_(radius),
      radius²_(radius_ * radius_) {}

template<typename Frame>
Position<Frame> const& Sphere<Frame>::centre() const {
  return centre_;
}

template<typename Frame>
Length const& Sphere<Frame>::radius() const {
  return radius_;
}

template<typename Frame>
Square<Length> const& Sphere<Frame>::radius²() const {
  return radius²_;
}

}  // namespace internal_sphere
}  // namespace geometry
}  // namespace principia
