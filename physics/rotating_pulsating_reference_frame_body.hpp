#pragma once

#include "physics/rotating_pulsating_reference_frame.hpp"

namespace principia {
namespace physics {
namespace _rotating_pulsating_reference_frame {
namespace internal {

using namespace principia::geometry::_homothecy;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

template<typename InertialFrame, typename ThisFrame>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
    RotatingPulsatingReferenceFrame(
        not_null<Ephemeris<InertialFrame> const*> ephemeris,
        not_null<MassiveBody const*> primary,
        not_null<MassiveBody const*> secondary)
    : ephemeris_(ephemeris),
      primary_(primary),
      secondary_(secondary),
      rotating_frame_(ephemeris_, primary_, secondary_) {}

template<typename InertialFrame, typename ThisFrame>
not_null<MassiveBody const*>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::primary() const {
  return primary_;
}

template<typename InertialFrame, typename ThisFrame>
not_null<MassiveBody const*>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::secondary() const {
  return secondary_;
}

template<typename InertialFrame, typename ThisFrame>
Instant RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::t_min()
    const {
  return rotating_frame_.t_min();
}

template<typename InertialFrame, typename ThisFrame>
Instant RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::t_max()
    const {
  return rotating_frame_.t_max();
}

template<typename InertialFrame, typename ThisFrame>
SimilarMotion<InertialFrame, ThisFrame> RotatingPulsatingReferenceFrame<
    InertialFrame,
    ThisFrame>::ToThisFrameAtTimeSimilarly(Instant const& t) const {
  return ToRotatingFrame(t).Inverse() *
         rotating_frame_.ToThisFrameAtTimeSimilarly(t);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame> RotatingPulsatingReferenceFrame<
    InertialFrame,
    ThisFrame>::GeometricAcceleration(Instant const& t,
                                      DegreesOfFreedom<ThisFrame> const&
                                          degrees_of_freedom) const {
  Displacement<InertialFrame> const u =
      ephemeris_->trajectory(primary_).EvaluatePosition(t) -
      ephemeris_->trajectory(secondary_).EvaluatePosition(t);
  Velocity<InertialFrame> const v =
      ephemeris_->trajectory(primary_).EvaluateVelocity(t) -
      ephemeris_->trajectory(secondary_).EvaluateVelocity(t);
  Vector<Acceleration, InertialFrame> const γ =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(primary_, t) -
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_, t);
  Length const r = u.Norm();
  Speed const ṙ = InnerProduct(u, v) / r;
  Speed const r̈ = v.Norm²() / r + InnerProduct(u, γ) / r - Pow<2>(ṙ) / r;
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame(t);
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  Vector<Acceleration, RotatingFrame> q̈ᴿ =
      rotating_frame_.GeometricAcceleration(
          t, to_rotating_frame(degrees_of_freedom));
  Displacement<ThisFrame> qᴾ =
      degrees_of_freedom.position() - ThisFrame::origin;
  Velocity<ThisFrame> q̇ᴾ = degrees_of_freedom.velocity();
  // (4.3)
  return -r̈ / r * qᴾ - 2 * ṙ / r * q̇ᴾ + from_rotating_frame.conformal_map()(q̈ᴿ);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
    RotationFreeGeometricAccelerationAtRest(
        Instant const& t,
        Position<ThisFrame> const& position) const {
  return Vector<Acceleration, ThisFrame>();
}

template<typename InertialFrame, typename ThisFrame>
SpecificEnergy
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::GeometricPotential(
    Instant const& t,
    Position<ThisFrame> const& position) const {
  return SpecificEnergy();
}

template<typename InertialFrame, typename ThisFrame>
void RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::WriteToMessage(
    not_null<serialization::ReferenceFrame*> message) const {
  LOG(FATAL) << "NOT IMPLEMENTED";
}

template<typename InertialFrame, typename ThisFrame>
auto RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::ToRotatingFrame(
    Instant const& t) const -> SimilarMotion<ThisFrame, RotatingFrame> {
  Length const r = (ephemeris_->trajectory(primary_).EvaluatePosition(t) -
                    ephemeris_->trajectory(secondary_).EvaluatePosition(t))
                       .Norm();
  Speed const ṙ = (ephemeris_->trajectory(primary_).EvaluateVelocity(t) -
                   ephemeris_->trajectory(secondary_).EvaluateVelocity(t))
                      .Norm();
  return SimilarMotion<ThisFrame, RotatingFrame>::DilatationAboutOrigin(
      Homothecy<double, ThisFrame, RotatingFrame>(r / 1 * Metre), ṙ / r)
}

}  // namespace internal
}  // namespace _rotating_pulsating_reference_frame
}  // namespace physics
}  // namespace principia