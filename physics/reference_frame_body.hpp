#pragma once

#include "physics/reference_frame.hpp"

#include "physics/rigid_reference_frame.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _reference_frame {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

template<typename InertialFrame, typename ThisFrame>
SimilarMotion<InertialFrame, ThisFrame>
ReferenceFrame<InertialFrame, ThisFrame>::ToThisFrameAtTimeSimilarly(
    Instant const& t) const {
  return FromThisFrameAtTimeSimilarly(t).Inverse();
}

template<typename InertialFrame, typename ThisFrame>
SimilarMotion<ThisFrame, InertialFrame>
ReferenceFrame<InertialFrame, ThisFrame>::FromThisFrameAtTimeSimilarly(
    Instant const& t) const {
  return ToThisFrameAtTimeSimilarly(t).Inverse();
}

template<typename InertialFrame, typename ThisFrame>
Rotation<Frenet<ThisFrame>, ThisFrame>
ReferenceFrame<InertialFrame, ThisFrame>::FrenetFrame(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Velocity<ThisFrame> const& velocity = degrees_of_freedom.velocity();
  Vector<Acceleration, ThisFrame> const acceleration =
      GeometricAcceleration(t, degrees_of_freedom);
  Vector<Acceleration, ThisFrame> const normal_acceleration =
      acceleration.OrthogonalizationAgainst(velocity);
  Vector<double, ThisFrame> tangent = Normalize(velocity);
  Vector<double, ThisFrame> normal = Normalize(normal_acceleration);
  Bivector<double, ThisFrame> binormal = Wedge(tangent, normal);
  // Maps |tangent| to {1, 0, 0}, |normal| to {0, 1, 0}, and |binormal| to
  // {0, 0, 1}.
  return Rotation<Frenet<ThisFrame>, ThisFrame>(tangent, normal, binormal);
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<ReferenceFrame<InertialFrame, ThisFrame>>>
ReferenceFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    serialization::ReferenceFrame const& message,
    not_null<Ephemeris<InertialFrame> const*> const ephemeris) {
  if (message.HasExtension(
          serialization::RotatingPulsatingReferenceFrame::extension)) {
    return static_cast<not_null<std::unique_ptr<ReferenceFrame>>>(
        RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(ephemeris,
                            message.GetExtension(
                                serialization::RotatingPulsatingReferenceFrame::
                                    extension)));
  } else {
    return static_cast<not_null<std::unique_ptr<ReferenceFrame>>>(
        RigidReferenceFrame<InertialFrame, ThisFrame>::ReadFromMessage(
            message, ephemeris));
  }
}

}  // namespace internal
}  // namespace _reference_frame
}  // namespace physics
}  // namespace principia
