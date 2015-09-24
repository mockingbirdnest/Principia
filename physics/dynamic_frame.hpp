#pragma once

#include <functional>

#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::AffineMap;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;
using quantities::Length;

namespace physics {

template<typename FromFrame, typename ToFrame>
class InstantaneousFrameMap {
 public:
  explicit InstantaneousFrameMap(
      AffineMap<FromFrame, ToFrame, Length, Rotation> const& map);

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;
  Rotation<FromFrame, ToFrame> const& linear_map() const;
  AffineMap<FromFrame, ToFrame, Length, Rotation> const& affine_map() const;

 private:
  AffineMap<FromFrame, ToFrame, Length, Rotation> map_;
};

// The definition of a reference frame |Frame| in arbitrary motion with respect
// to |BaseFrame|.
template<typename BaseFrame, typename Frame>
class DynamicFrame {
 public:
  virtual InstantaneousFrameMap<BaseFrame, Frame> ToFrameAtTime(
      Instant const& t) = 0;
  virtual InstantaneousFrameMap<Frame, BaseFrame> FromFrameAtTime(
      Instant const& t) = 0;

  // The acceleration due to the non-inertial motion of |Frame| and gravity.
  virtual Vector<Acceleration, Frame> GeometricAcceleration(
      Instant const& t,
      Position const& q) = 0;
};

}  // namespace physics
}  // namespace principia
