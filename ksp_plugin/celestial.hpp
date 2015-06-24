#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/mobile_interface.hpp"
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::not_null;
using physics::Body;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using quantities::GravitationalParameter;

namespace ksp_plugin {

// Represents a KSP |CelestialBody|.
class Celestial {
 public:
  explicit Celestial(not_null<MassiveBody const*> body);
  Celestial(Celestial const&) = delete;
  Celestial(Celestial&&) = delete;
  ~Celestial() = default;

  // True if, and only if, |trajectory_| is not null.
  bool is_initialized() const;
  void set_trajectory(
      not_null<ContinuousTrajectory<Barycentric> const*> const trajectory);
  ContinuousTrajectory<Barycentric> const& trajectory() const;
  not_null<ContinuousTrajectory<Barycentric>::Hint*> current_time_hint() const;
  DegreesOfFreedom<Barycentric> current_degrees_of_freedom(
      Instant const& current_time) const;
  Position<Barycentric> current_position(Instant const& current_time) const;
  Velocity<Barycentric> current_velocity(Instant const& current_time) const;

  MassiveBody const& body() const;
  bool has_parent() const;
  Celestial const* parent() const;  // Null for the Sun.
  void set_parent(not_null<Celestial const*> const parent);

 private:
  not_null<MassiveBody const*> body_;
  // The parent body for the 2-body approximation. Not owning, must only
  // be null for the sun.
  Celestial const* parent_ = nullptr;
  ContinuousTrajectory<Barycentric> const* trajectory_ = nullptr;
  not_null<
      std::unique_ptr<
          ContinuousTrajectory<Barycentric>::Hint>> current_time_hint_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/celestial_body.hpp"
