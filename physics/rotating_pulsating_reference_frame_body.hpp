#pragma once

#include "physics/rotating_pulsating_reference_frame.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "geometry/homothecy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _rotating_pulsating_reference_frame {
namespace internal {

using namespace principia::geometry::_homothecy;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

template<typename InertialFrame, typename ThisFrame>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
    RotatingPulsatingReferenceFrame(
        not_null<Ephemeris<InertialFrame> const*> const ephemeris,
        not_null<MassiveBody const*> const primary,
        not_null<MassiveBody const*> const secondary)
    : RotatingPulsatingReferenceFrame(
          ephemeris,
          std::vector<not_null<MassiveBody const*>>{primary},
          std::vector<not_null<MassiveBody const*>>{secondary}) {}

template<typename InertialFrame, typename ThisFrame>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
    RotatingPulsatingReferenceFrame(
        not_null<Ephemeris<InertialFrame> const*> const ephemeris,
        std::vector<not_null<MassiveBody const*>> primaries,
        std::vector<not_null<MassiveBody const*>> secondaries)
    : ephemeris_(ephemeris),
      primaries_(std::move(primaries)),
      secondaries_(std::move(secondaries)),
      rotating_frame_(ephemeris_, primaries_, secondaries_) {
  absl::btree_set<not_null<MassiveBody const*>> primary_set(primaries_.begin(),
                                                            primaries_.end());
  absl::btree_set<not_null<MassiveBody const*>> secondary_set(
      secondaries_.begin(), secondaries_.end());
  absl::btree_set<not_null<MassiveBody const*>> intersection;
  std::set_intersection(primary_set.begin(),
                        primary_set.end(),
                        secondary_set.begin(),
                        secondary_set.end(),
                        std::inserter(intersection, intersection.begin()));
  auto const names = [](auto const& bodies) {
    return absl::StrJoin(
        bodies,
        ",",
        [](std::string* const out, not_null<MassiveBody const*> const body) {
          out->append(body->name());
        });
  };
  CHECK_GE(primaries_.size(), 1) << names(primaries_);
  CHECK_EQ(primary_set.size(), primaries_.size()) << names(primaries_);
  CHECK_GE(secondaries_.size(), 1) << names(secondaries_);
  CHECK_EQ(secondary_set.size(), secondaries_.size()) << names(secondaries_);
  CHECK_EQ(intersection.size(), 0) << names(intersection);
}

template<typename InertialFrame, typename ThisFrame>
std::vector<not_null<MassiveBody const*>> const&
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::primaries() const {
  return primaries_;
}

template<typename InertialFrame, typename ThisFrame>
std::vector<not_null<MassiveBody const*>> const&
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::secondaries() const {
  return secondaries_;
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
  return ToRotatingFrame(r_derivatives<1>(t)).Inverse() *
         rotating_frame_.ToThisFrameAtTimeSimilarly(t);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame> RotatingPulsatingReferenceFrame<
    InertialFrame,
    ThisFrame>::GeometricAcceleration(Instant const& t,
                                      DegreesOfFreedom<ThisFrame> const&
                                          degrees_of_freedom) const {
  auto const [r, ṙ, r̈] = r_derivatives<2>(t);
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame({r, ṙ});
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  Vector<Acceleration, RotatingFrame> const q̈ᴿ =
      rotating_frame_.GeometricAcceleration(
          t, to_rotating_frame(degrees_of_freedom));
  Displacement<ThisFrame> const qᴾ =
      degrees_of_freedom.position() - ThisFrame::origin;
  Velocity<ThisFrame> const q̇ᴾ = degrees_of_freedom.velocity();
  // See equation (4.3) in Rotating Pulsating.pdf.
  return -r̈ / r * qᴾ - 2 * ṙ / r * q̇ᴾ + from_rotating_frame.conformal_map()(q̈ᴿ);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::
    RotationFreeGeometricAccelerationAtRest(
        Instant const& t,
        Position<ThisFrame> const& position) const {
  auto const [r, ṙ, r̈] = r_derivatives<2>(t);
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame({r, ṙ});
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  Vector<Acceleration, RotatingFrame> const dVᴿⳆdqᴿ =
      rotating_frame_.RotationFreeGeometricAccelerationAtRest(
          t, to_rotating_frame.similarity()(position));
  Displacement<ThisFrame> const qᴾ = position - ThisFrame::origin;
  // See equations (4.3) and (4.4) in Rotating Pulsating.pdf.
  return -r̈ / r * qᴾ + from_rotating_frame.conformal_map()(dVᴿⳆdqᴿ);
}

template<typename InertialFrame, typename ThisFrame>
SpecificEnergy
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::GeometricPotential(
    Instant const& t,
    Position<ThisFrame> const& position) const {
  auto const [r, ṙ, r̈] = r_derivatives<2>(t);
  SimilarMotion<ThisFrame, RotatingFrame> const to_rotating_frame =
      ToRotatingFrame({r, ṙ});
  SimilarMotion<RotatingFrame, ThisFrame> const from_rotating_frame =
      to_rotating_frame.Inverse();
  SpecificEnergy const Vᴿ = rotating_frame_.GeometricPotential(
      t, to_rotating_frame.similarity()(position));
  Displacement<ThisFrame> const qᴾ = position - ThisFrame::origin;
  // See Vᴾ in equation (4.4) in Rotating Pulsating.pdf.
  return r̈ * qᴾ.Norm²() / (2 * r) + Vᴿ / Pow<2>(r / (1 * Metre));
}

template<typename InertialFrame, typename ThisFrame>
void RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::WriteToMessage(
    not_null<serialization::ReferenceFrame*> message) const {
  auto* const extension = message->MutableExtension(
      serialization::RotatingPulsatingReferenceFrame::extension);
  for (not_null const primary : primaries_) {
    extension->add_primary(
        ephemeris_->serialization_index_for_body(primary));
  }
  for (not_null const secondary : secondaries_) {
    extension->add_secondary(
        ephemeris_->serialization_index_for_body(secondary));
  }
}

template<typename InertialFrame, typename ThisFrame>
not_null<
    std::unique_ptr<RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>>>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::RotatingPulsatingReferenceFrame const& message) {
  std::vector<not_null<MassiveBody const*>> primaries;
  primaries.reserve(message.primary().size());
  for (int const primary : message.primary()) {
    primaries.push_back(ephemeris->body_for_serialization_index(primary));
  }
  std::vector<not_null<MassiveBody const*>> secondaries;
  secondaries.reserve(message.secondary().size());
  for (int const secondary : message.secondary()) {
    secondaries.push_back(ephemeris->body_for_serialization_index(secondary));
  }
  return std::make_unique<RotatingPulsatingReferenceFrame>(
      ephemeris, std::move(primaries), std::move(secondaries));
}

template<typename InertialFrame, typename ThisFrame>
template<int degree>
Derivatives<Length, Instant, degree + 1>
RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::r_derivatives(
    Instant const& t) const {
  Displacement<InertialFrame> const u =
      rotating_frame_.template PrimaryDerivative<0>(t) -
      rotating_frame_.template SecondaryDerivative<0>(t);
  Length const r = u.Norm();
  if constexpr (degree == 0) {
    return {r};
  } else {
    Velocity<InertialFrame> const v =
        rotating_frame_.template PrimaryDerivative<1>(t) -
        rotating_frame_.template SecondaryDerivative<1>(t);
    Speed const ṙ = InnerProduct(u, v) / r;
    if constexpr (degree == 1) {
      return {r, ṙ};
    } else {
      static_assert(degree == 2);
      Vector<Acceleration, InertialFrame> const γ =
          rotating_frame_.template PrimaryDerivative<2>(t) -
          rotating_frame_.template SecondaryDerivative<2>(t);
      Acceleration const r̈ =
          v.Norm²() / r + InnerProduct(u, γ) / r - Pow<2>(ṙ) / r;
      return {r, ṙ, r̈};
    }
  }
}

template<typename InertialFrame, typename ThisFrame>
auto RotatingPulsatingReferenceFrame<InertialFrame, ThisFrame>::ToRotatingFrame(
    Derivatives<Length, Instant, 2> const& r_derivatives_1) const
    -> SimilarMotion<ThisFrame, RotatingFrame> {
  auto const& [r, ṙ] = r_derivatives_1;
  return SimilarMotion<ThisFrame, RotatingFrame>::DilatationAboutOrigin(
      Homothecy<double, ThisFrame, RotatingFrame>(r / (1 * Metre)), ṙ / r);
}

}  // namespace internal
}  // namespace _rotating_pulsating_reference_frame
}  // namespace physics
}  // namespace principia
