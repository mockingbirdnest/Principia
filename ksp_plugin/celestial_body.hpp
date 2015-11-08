#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

inline Celestial::Celestial(not_null<MassiveBody const*> body)
    : body_(body),
      current_time_hint_(
          make_not_null_unique<ContinuousTrajectory<Barycentric>::Hint>()) {}

inline bool Celestial::is_initialized() const {
  return trajectory_ != nullptr;
}

inline void Celestial::set_trajectory(
      not_null<ContinuousTrajectory<Barycentric> const*> const trajectory) {
  CHECK(!is_initialized());
  trajectory_ = trajectory;
}

inline ContinuousTrajectory<Barycentric> const& Celestial::trajectory() const {
  CHECK(is_initialized());
  return *trajectory_;
}

inline not_null<ContinuousTrajectory<Barycentric>::Hint*>
Celestial::current_time_hint() const {
  CHECK(is_initialized());
  return current_time_hint_.get();
}

inline DegreesOfFreedom<Barycentric> Celestial::current_degrees_of_freedom(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluateDegreesOfFreedom(current_time,
                                               current_time_hint());
}

inline Position<Barycentric> Celestial::current_position(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluatePosition(current_time, current_time_hint());
}

inline Velocity<Barycentric> Celestial::current_velocity(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluateVelocity(current_time, current_time_hint());
}

inline not_null<MassiveBody const*> Celestial::body() const {
  return body_;
}

inline bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

inline Celestial const* Celestial::parent() const {
  return parent_;
}

inline void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

}  // namespace ksp_plugin
}  // namespace principia
