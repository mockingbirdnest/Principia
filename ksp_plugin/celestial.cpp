#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_celestial {

Celestial::Celestial(not_null<RotatingBody<Barycentric> const*> body)
    : body_(body),
      current_time_hint_(
          make_not_null_unique<ContinuousTrajectory<Barycentric>::Hint>()) {}

bool Celestial::is_initialized() const {
  return trajectory_ != nullptr;
}

void Celestial::set_trajectory(
    not_null<ContinuousTrajectory<Barycentric> const*> const trajectory) {
  CHECK(!is_initialized());
  trajectory_ = trajectory;
}

ContinuousTrajectory<Barycentric> const& Celestial::trajectory() const {
  CHECK(is_initialized());
  return *trajectory_;
}

not_null<ContinuousTrajectory<Barycentric>::Hint*>
Celestial::current_time_hint() const {
  CHECK(is_initialized());
  return current_time_hint_.get();
}

DegreesOfFreedom<Barycentric> Celestial::current_degrees_of_freedom(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluateDegreesOfFreedom(current_time,
                                               current_time_hint());
}

Position<Barycentric> Celestial::current_position(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluatePosition(current_time, current_time_hint());
}

Velocity<Barycentric> Celestial::current_velocity(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluateVelocity(current_time, current_time_hint());
}

not_null<RotatingBody<Barycentric> const*> Celestial::body() const {
  return body_;
}

bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

Celestial const* Celestial::parent() const {
  return parent_;
}

void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

}  // namespace internal_celestial
}  // namespace ksp_plugin
}  // namespace principia
