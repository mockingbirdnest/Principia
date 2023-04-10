#pragma once

#include "physics/reference_frame.hpp"

#include "physics/body_centred_body_direction_reference_frame.hpp"

namespace principia {
namespace physics {
namespace _rotating_pulsating_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_space;
using namespace principia::physics::_barycentric_rotating_reference_frame;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_tuples;

template<typename InertialFrame, typename ThisFrame>
class RotatingPulsatingReferenceFrame
    : public ReferenceFrame<InertialFrame, ThisFrame> {
  static_assert(ThisFrame::may_rotate);
  // TODO(egg): This frame may dilate, too.

 public:
  RotatingPulsatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> primary,
      not_null<MassiveBody const*> secondary);

  not_null<MassiveBody const*> primary() const;
  not_null<MassiveBody const*> secondary() const;

  Instant t_min() const override;
  Instant t_max() const override;

  SimilarMotion<InertialFrame, ThisFrame> ToThisFrameAtTimeSimilarly(
      Instant const& t) const override;

  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

  Vector<Acceleration, ThisFrame> RotationFreeGeometricAccelerationAtRest(
      Instant const& t,
      Position<ThisFrame> const& position) const override;

  SpecificEnergy GeometricPotential(
      Instant const& t,
      Position<ThisFrame> const& position) const override;

  void WriteToMessage(
      not_null<serialization::ReferenceFrame*> message) const override;

 private:
  using RotatingFrame =
      Frame<struct RotatingTag, Arbitrary, ThisFrame::handedness>;

  template<int degree>
  Derivatives<Length, Instant, degree + 1> r_taylor(Instant const& t) const;


  SimilarMotion<ThisFrame, RotatingFrame> ToRotatingFrame(
      Derivatives<Length, Instant, 2> const& r_taylor_1) const;

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  not_null<MassiveBody const*> const primary_;
  not_null<MassiveBody const*> const secondary_;
  BarycentricRotatingReferenceFrame<InertialFrame, RotatingFrame> const
      rotating_frame_;
};

}  // namespace internal

using internal::RotatingPulsatingReferenceFrame;

}  // namespace _rotating_pulsating_reference_frame
}  // namespace physics
}  // namespace principia

#include "physics/rotating_pulsating_reference_frame_body.hpp"
