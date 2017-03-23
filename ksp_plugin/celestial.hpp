
#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_celestial {

using base::not_null;
using base::make_not_null_unique;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::Body;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::RotatingBody;
using quantities::GravitationalParameter;

// Represents a KSP |CelestialBody|.
class Celestial final {
 public:
  explicit Celestial(not_null<RotatingBody<Barycentric> const*> body);
  Celestial(Celestial const&) = delete;
  Celestial(Celestial&&) = delete;

  // True if, and only if, |trajectory_| is not null.
  bool is_initialized() const;
  void set_trajectory(
      not_null<ContinuousTrajectory<Barycentric> const*> trajectory);
  ContinuousTrajectory<Barycentric> const& trajectory() const;
  not_null<ContinuousTrajectory<Barycentric>::Hint*> current_time_hint() const;
  DegreesOfFreedom<Barycentric> current_degrees_of_freedom(
      Instant const& current_time) const;
  Position<Barycentric> current_position(Instant const& current_time) const;
  Velocity<Barycentric> current_velocity(Instant const& current_time) const;

  not_null<RotatingBody<Barycentric> const*> body() const;
  bool has_parent() const;
  Celestial const* parent() const;  // Null for the Sun.
  void set_parent(not_null<Celestial const*> parent);

 private:
  not_null<RotatingBody<Barycentric> const*> body_;
  // The parent body for the 2-body approximation. Not owning, must only
  // be null for the sun.
  Celestial const* parent_ = nullptr;
  ContinuousTrajectory<Barycentric> const* trajectory_ = nullptr;
  not_null<
      std::unique_ptr<
          ContinuousTrajectory<Barycentric>::Hint>> current_time_hint_;
};

}  // namespace internal_celestial

using internal_celestial::Celestial;

}  // namespace ksp_plugin
}  // namespace principia
