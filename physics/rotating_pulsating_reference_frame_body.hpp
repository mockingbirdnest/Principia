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
  return ToRotatingFrame(r_taylor<1>(t)).Inverse() *
         rotating_frame_.ToThisFrameAtTimeSimilarly(t);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame> RotatingPulsatingReferenceFrame<
    InertialFrame,
    ThisFrame>::GeometricAcceleration(Instant const& t,
                                      DegreesOfFreedom<ThisFrame> const&
                                          degrees_of_freedom) const {
  auto const [r, ṙ, r̈] = r_taylor<2>(t);
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame({r, ṙ});
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  Vector<Acceleration, RotatingFrame> q̈ᴿ =
      rotating_frame_.GeometricAcceleration(
          t, to_rotating_frame(degrees_of_freedom));
  Displacement<ThisFrame> qᴾ =
      degrees_of_freedom.position() - ThisFrame::origin;
  Velocity<ThisFrame> q̇ᴾ = degrees_of_freedom.velocity();
  // See equation (4.3).
  return -r̈ / r * qᴾ - 2 * ṙ / r * q̇ᴾ + from_rotating_frame.conformal_map()(q̈ᴿ);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
    RotationFreeGeometricAccelerationAtRest(
        Instant const& t,
        Position<ThisFrame> const& position) const {
  auto const [r, ṙ, r̈] = r_taylor<2>(t);
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame({r, ṙ});
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  Vector<Acceleration, RotatingFrame> dVᴿⳆdqᴿ =
      rotating_frame_.RotationFreeGeometricAccelerationAtRest(
          t, to_rotating_frame.similarity()(position));
  Displacement<ThisFrame> qᴾ = position - ThisFrame::origin;
  // See equations (4.3) and (4.4).
  return -r̈ / r * qᴾ + from_rotating_frame.conformal_map()(dVᴿⳆdqᴿ);
}

template<typename InertialFrame, typename ThisFrame>
SpecificEnergy
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::GeometricPotential(
    Instant const& t,
    Position<ThisFrame> const& position) const {
  auto const [r, ṙ, r̈] = r_taylor<2>(t);
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame({r, ṙ});
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  SpecificEnergy Vᴿ = rotating_frame_.GeometricPotential(
      t, to_rotating_frame.similarity()(position));
  Displacement<ThisFrame> qᴾ = position - ThisFrame::origin;
  // See Vᴾ in equation (4.4).
  return r̈ * qᴾ.Norm²() / (2 * r) + Vᴿ / Pow<2>(r / (1 * Metre));
}

template<typename InertialFrame, typename ThisFrame>
void RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::WriteToMessage(
    not_null<serialization::ReferenceFrame*> message) const {
  LOG(FATAL) << "NOT IMPLEMENTED";
}

template<typename InertialFrame, typename ThisFrame>
template<int degree>
Derivatives<Length, Instant, degree + 1>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::r_taylor(
    Instant const& t) const {
  Displacement<InertialFrame> const u =
      ephemeris_->trajectory(primary_)->EvaluatePosition(t) -
      ephemeris_->trajectory(secondary_)->EvaluatePosition(t);
  Length const r = u.Norm();
  if constexpr (degree == 0) {
    return {r};
  } else {
    Velocity<InertialFrame> const v =
        ephemeris_->trajectory(primary_)->EvaluateVelocity(t) -
        ephemeris_->trajectory(secondary_)->EvaluateVelocity(t);
    Speed const ṙ = InnerProduct(u, v) / r;
    if constexpr (degree == 1) {
      return {r, ṙ};
    } else {
      Vector<Acceleration, InertialFrame> const γ =
          ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(primary_,
                                                                    t) -
          ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_,
                                                                    t);
      Acceleration const r̈ =
          v.Norm²() / r + InnerProduct(u, γ) / r - Pow<2>(ṙ) / r;
      static_assert(degree == 2);
      return {r, ṙ, r̈};
    }
  }
}

template<typename InertialFrame, typename ThisFrame>
auto RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::ToRotatingFrame(
    Derivatives<Length, Instant, 2> const& r_taylor_1) const
    -> SimilarMotion<ThisFrame, RotatingFrame> {
  auto const& [r, ṙ] = r_taylor_1;
  return SimilarMotion<ThisFrame, RotatingFrame>::DilatationAboutOrigin(
      Homothecy<double, ThisFrame, RotatingFrame>(r / (1 * Metre)), ṙ / r);
}

}  // namespace internal
}  // namespace _rotating_pulsating_reference_frame
}  // namespace physics
}  // namespace principia