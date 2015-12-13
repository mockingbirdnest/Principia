#pragma once

#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/body_centered_non_rotating_dynamic_frame.hpp"
#include "physics/dynamic_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Bivector;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::R3x3Matrix;
using geometry::Wedge;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Second;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
Rotation<Frenet<ThisFrame>, ThisFrame>
DynamicFrame<InertialFrame, ThisFrame>::FrenetFrame(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Velocity<ThisFrame> const& velocity = degrees_of_freedom.velocity();
  Vector<Acceleration, ThisFrame> const acceleration =
      GeometricAcceleration(t, degrees_of_freedom);
  Vector<Acceleration, ThisFrame> normal_acceleration = acceleration;
  velocity.template Orthogonalize<Acceleration>(&normal_acceleration);
  Vector<double, ThisFrame> tangent = Normalize(velocity);
  Vector<double, ThisFrame> normal = Normalize(normal_acceleration);
  Bivector<double, ThisFrame> binormal = Wedge(tangent, normal);
  // Maps |tangent| to {1, 0, 0}, |normal| to {0, 1, 0}, and |binormal| to
  // {0, 0, 1}.
  return Rotation<Frenet<ThisFrame>, ThisFrame>(
      R3x3Matrix(tangent.coordinates(),
                 normal.coordinates(),
                 binormal.coordinates()).Transpose());
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<DynamicFrame<InertialFrame, ThisFrame>>>
DynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::DynamicFrame const& message) {
  int extensions_found = 0;
  if (message.HasExtension(serialization::BarycentricRotatingDynamicFrame::
                               barycentric_rotating_dynamic_frame)) {
    ++extensions_found;
    return BarycentricRotatingDynamicFrame::ReadFromMessage(
               ephemeris,
               message.GetExtension(
                  serialization::BarycentricRotatingDynamicFrame::
                      barycentric_rotating_dynamic_frame));
  }
  if (message.HasExtension(serialization::BodyCentredNonRotatingDynamicFrame::
                               body_centred_non_rotating_dynamic_frame)) {
    ++extensions_found;
    return BodyCentredNonRotatingDynamicFrame::ReadFromMessage(
               ephemeris,
               message.GetExtension(
                   serialization::BodyCentredNonRotatingDynamicFrame::
                       body_centred_non_rotating_dynamic_frame));
  }
  CHECK_EQ(1, extensions_found) << message.DebugString();
}

}  // namespace physics
}  // namespace principia
