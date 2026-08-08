#pragma once

#include "physics/reference_frame.hpp"

#include <memory>

#include "geometry/r3x3_matrix.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _reference_frame {
namespace internal {

using namespace principia::geometry::_r3x3_matrix;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::physics::_rotating_pulsating_reference_frame;
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
  Vector<double, ThisFrame> const tangent = Normalize(velocity);
  Vector<double, ThisFrame> const normal = Normalize(normal_acceleration);
  Bivector<double, ThisFrame> const binormal = Wedge(tangent, normal);
  // Maps `tangent` to {1, 0, 0}, `normal` to {0, 1, 0}, and `binormal` to
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

template<typename InertialFrame, typename ThisFrame>
bool operator==(ReferenceFrame<InertialFrame, ThisFrame> const& left,
                ReferenceFrame<InertialFrame, ThisFrame> const& right) {
  // We exceptionally use RTTI here because `operator==` is hopelessly broken
  // in C++ in the presence of inheritance.
  if (typeid(left) != typeid(right)) {
    return false;
  }
  {
    using RotatingPulsating =
        RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>;
    if (typeid(left) == typeid(RotatingPulsating)) {
      return static_cast<RotatingPulsating const&>(left) ==
             static_cast<RotatingPulsating const&>(right);
    }
  }
  {
    // The type must be descended from `RigidReferenceFrame`, but we cannot
    // check this using `typeid`.
    // TODO(phl): Improve once we have reflection in C++26.
    using Rigid = RigidReferenceFrame<InertialFrame, ThisFrame>;
    auto const* left_rigid = dynamic_cast<Rigid const*>(&left);
    auto const* right_rigid = dynamic_cast<Rigid const*>(&right);
    CHECK_NE(nullptr, left_rigid)
        << "The type of left is " << typeid(left).name();
    CHECK_NE(nullptr, right_rigid)
        << "The type of right is " << typeid(right).name();
    return *left_rigid == *right_rigid;
  }
}

}  // namespace internal
}  // namespace _reference_frame
}  // namespace physics
}  // namespace principia
