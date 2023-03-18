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
namespace _celestial {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_named_quantities;
using namespace principia::physics::_body;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_named_quantities;

// Represents a KSP |CelestialBody|.
class Celestial {
 public:
  explicit Celestial(not_null<RotatingBody<Barycentric> const*> body);
  Celestial(Celestial const&) = delete;
  Celestial(Celestial&&) = delete;
  virtual ~Celestial() = default;

  // True if, and only if, |trajectory_| is not null.
  bool is_initialized() const;
  void set_trajectory(
      not_null<ContinuousTrajectory<Barycentric> const*> trajectory);
  ContinuousTrajectory<Barycentric> const& trajectory() const;
  virtual DegreesOfFreedom<Barycentric> current_degrees_of_freedom(
      Instant const& current_time) const;
  virtual Position<Barycentric> current_position(
      Instant const& current_time) const;
  virtual Velocity<Barycentric> current_velocity(
      Instant const& current_time) const;

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
};

}  // namespace internal

using internal::Celestial;

}  // namespace _celestial
}  // namespace ksp_plugin
}  // namespace principia

namespace principia::ksp_plugin {
using namespace principia::ksp_plugin::_celestial;
}  // namespace principia::ksp_plugin
