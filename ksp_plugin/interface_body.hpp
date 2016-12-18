
#pragma once

#include "ksp_plugin/interface.hpp"

#include <cmath>
#include <limits>

namespace principia {
namespace interface {

using geometry::Position;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using physics::Ephemeris;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

inline bool NaNIndependentEq(double const left, double const right) {
  return (left == right) || (std::isnan(left) && std::isnan(right));
}

template<typename Container>
TypedIterator<Container>::TypedIterator(Container container)
    : container_(std::move(container)),
      iterator_(container_.begin()) {}

template<typename Container>
template<typename Interchange>
Interchange TypedIterator<Container>::Get(
    std::function<Interchange(typename Container::value_type const&)> const&
        convert) const {
  CHECK(iterator_ != container_.end());
  return convert(*iterator_);
}

template<typename Container>
bool TypedIterator<Container>::AtEnd() const {
  return iterator_ == container_.end();
}

template<typename Container>
void TypedIterator<Container>::Increment() {
  ++iterator_;
}

template<typename Container>
int TypedIterator<Container>::Size() const {
  return container_.size();
}

inline TypedIterator<DiscreteTrajectory<World>>::TypedIterator(
    not_null<std::unique_ptr<DiscreteTrajectory<World>>> trajectory,
    not_null<Plugin const*> const plugin)
    : trajectory_(std::move(trajectory)),
      iterator_(trajectory_->Begin()),
      plugin_(plugin) {
  CHECK(trajectory_->is_root());
}

template<typename Interchange>
Interchange TypedIterator<DiscreteTrajectory<World>>::Get(
    std::function<Interchange(
        DiscreteTrajectory<World>::Iterator const&)> const& convert) const {
  CHECK(iterator_ != trajectory_->End());
  return convert(iterator_);
}

inline bool TypedIterator<DiscreteTrajectory<World>>::AtEnd() const {
  return iterator_ == trajectory_->End();
}

inline void TypedIterator<DiscreteTrajectory<World>>::Increment() {
  ++iterator_;
}

inline int TypedIterator<DiscreteTrajectory<World>>::Size() const {
  return trajectory_->Size();
}

inline not_null<Plugin const*> TypedIterator<
    DiscreteTrajectory<World>>::plugin() const {
  return plugin_;
}


template<typename T>
std::unique_ptr<T> TakeOwnership(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

template<typename T>
std::unique_ptr<T[]> TakeOwnershipArray(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T[]> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

inline bool operator==(AdaptiveStepParameters const& left,
                       AdaptiveStepParameters const& right) {
  return left.max_steps == right.max_steps &&
         NaNIndependentEq(left.length_integration_tolerance,
                          right.length_integration_tolerance) &&
         NaNIndependentEq(left.speed_integration_tolerance,
                          right.speed_integration_tolerance);
}

inline bool operator==(Burn const& left, Burn const& right) {
  return NaNIndependentEq(left.thrust_in_kilonewtons,
                          right.thrust_in_kilonewtons) &&
         NaNIndependentEq(left.specific_impulse_in_seconds_g0,
                          right.specific_impulse_in_seconds_g0) &&
         left.frame == right.frame &&
         NaNIndependentEq(left.initial_time, right.initial_time) &&
         left.delta_v == right.delta_v;
}

inline bool operator==(NavigationFrameParameters const& left,
                       NavigationFrameParameters const& right) {
  return left.extension == right.extension &&
         left.centre_index == right.centre_index &&
         left.primary_index == right.primary_index &&
         left.secondary_index == right.secondary_index;
}

inline bool operator==(NavigationManoeuvre const& left,
                       NavigationManoeuvre const& right) {
  return left.burn == right.burn &&
         NaNIndependentEq(left.initial_mass_in_tonnes,
                          right.initial_mass_in_tonnes) &&
         NaNIndependentEq(left.final_mass_in_tonnes,
                          right.final_mass_in_tonnes) &&
         NaNIndependentEq(left.mass_flow, right.mass_flow) &&
         NaNIndependentEq(left.duration, right.duration) &&
         NaNIndependentEq(left.final_time, right.final_time) &&
         NaNIndependentEq(left.time_of_half_delta_v,
                          right.time_of_half_delta_v) &&
         NaNIndependentEq(left.time_to_half_delta_v,
                          right.time_to_half_delta_v) &&
         left.inertial_direction == right.inertial_direction &&
         left.binormal == right.binormal &&
         left.normal == right.normal &&
         left.tangent == right.tangent;
}

inline bool operator==(QP const& left, QP const& right) {
  return left.q == right.q && left.p == right.p;
}

inline bool operator==(WXYZ const& left, WXYZ const& right) {
  return NaNIndependentEq(left.w, right.w) &&
         NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y) &&
         NaNIndependentEq(left.z, right.z);
}

inline bool operator==(XYZ const& left, XYZ const& right) {
  return NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y) &&
         NaNIndependentEq(left.z, right.z);
}

inline physics::Ephemeris<Barycentric>::AdaptiveStepParameters
FromAdaptiveStepParameters(
    AdaptiveStepParameters const& adaptive_step_parameters) {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
      adaptive_step_parameters.max_steps,
      adaptive_step_parameters.length_integration_tolerance * Metre,
      adaptive_step_parameters.speed_integration_tolerance * (Metre / Second));
}

inline physics::KeplerianElements<Barycentric> FromKeplerianElements(
    KeplerianElements const& keplerian_elements) {
  physics::KeplerianElements<Barycentric> barycentric_keplerian_elements;
  barycentric_keplerian_elements.eccentricity = keplerian_elements.eccentricity;
  barycentric_keplerian_elements.semimajor_axis =
      std::isnan(keplerian_elements.semimajor_axis)
          ? std::experimental::nullopt
          : std::experimental::make_optional(keplerian_elements.semimajor_axis *
                                             Metre);
  barycentric_keplerian_elements.mean_motion =
      std::isnan(keplerian_elements.mean_motion)
          ? std::experimental::nullopt
          : std::experimental::make_optional(keplerian_elements.mean_motion *
                                             Radian / Second);
  barycentric_keplerian_elements.inclination =
      keplerian_elements.inclination_in_degrees * Degree;
  barycentric_keplerian_elements.longitude_of_ascending_node =
      keplerian_elements.longitude_of_ascending_node_in_degrees * Degree;
  barycentric_keplerian_elements.argument_of_periapsis =
      keplerian_elements.argument_of_periapsis_in_degrees * Degree;
  barycentric_keplerian_elements.mean_anomaly =
      keplerian_elements.mean_anomaly * Radian;
  return barycentric_keplerian_elements;
}

inline geometry::R3Element<double> FromXYZ(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

inline AdaptiveStepParameters ToAdaptiveStepParameters(
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  return {adaptive_step_parameters.max_steps(),
          adaptive_step_parameters.length_integration_tolerance() / Metre,
          adaptive_step_parameters.speed_integration_tolerance() /
              (Metre / Second)};
}

inline KeplerianElements ToKeplerianElements(
    physics::KeplerianElements<Barycentric> const& keplerian_elements) {
  return {keplerian_elements.eccentricity,
          keplerian_elements.semimajor_axis
              ? *keplerian_elements.semimajor_axis / Metre
              : std::numeric_limits<double>::quiet_NaN(),
          keplerian_elements.mean_motion
              ? *keplerian_elements.mean_motion / (Radian / Second)
              : std::numeric_limits<double>::quiet_NaN(),
          keplerian_elements.inclination / Degree,
          keplerian_elements.longitude_of_ascending_node / Degree,
          keplerian_elements.argument_of_periapsis / Degree,
          keplerian_elements.mean_anomaly / Radian};
}

inline WXYZ ToWXYZ(geometry::Quaternion const& quaternion) {
  return {quaternion.real_part(),
          quaternion.imaginary_part().x,
          quaternion.imaginary_part().y,
          quaternion.imaginary_part().z};
}

inline XYZ ToXYZ(geometry::R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
}

inline Instant FromGameTime(Plugin const& plugin,
                            double const t) {
  return plugin.GameEpoch() + t * Second;
}

inline double ToGameTime(Plugin const& plugin,
                         Instant const& t) {
  return (t - plugin.GameEpoch()) / Second;
}

inline not_null<Vessel*> GetVessel(Plugin const& plugin,
                                   char const* const vessel_guid) {
  CHECK(plugin.HasVessel(vessel_guid)) << vessel_guid;
  return plugin.GetVessel(vessel_guid);
}

inline not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const& plugin,
    NavigationFrameParameters const& parameters) {
  switch (parameters.extension) {
    case serialization::BarycentricRotatingDynamicFrame::
        kExtensionFieldNumber:
      return plugin.NewBarycentricRotatingNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredBodyDirectionDynamicFrame::
        kExtensionFieldNumber:
      return plugin.NewBodyCentredBodyDirectionNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredNonRotatingDynamicFrame::
        kExtensionFieldNumber:
      return plugin.NewBodyCentredNonRotatingNavigationFrame(
          parameters.centre_index);
    case serialization::BodySurfaceDynamicFrame::kExtensionFieldNumber:
      return plugin.NewBodySurfaceNavigationFrame(parameters.centre_index);
    default:
      LOG(FATAL) << "Unexpected extension " << parameters.extension;
      base::noreturn();
  }
}

}  // namespace interface
}  // namespace principia
