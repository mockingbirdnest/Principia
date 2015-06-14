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
  ContinuousTrajectory<Barycentric> const& trajectory() const;
  not_null<ContinuousTrajectory<Barycentric>::Hint*> current_time_hint();

  MassiveBody const& body() const;
  bool has_parent() const;
  Celestial const* parent() const;  // Null for the Sun.
  void set_parent(not_null<Celestial const*> const parent);

  // The celestial must satisfy |is_initialized()|.
  void WriteToMessage(not_null<serialization::Celestial*> const message) const;
  // NOTE(egg): This should return a |not_null|, but we can't do that until
  // |not_null<std::unique_ptr<T>>| is convertible to |std::unique_ptr<T>|, and
  // that requires a VS 2015 feature (rvalue references for |*this|).
  static std::unique_ptr<Celestial> ReadFromMessage(
      serialization::Celestial const& message);

 private:
  not_null<MassiveBody const*> body_;
  // The parent body for the 2-body approximation. Not owning, must only
  // be null for the sun.
  Celestial const* parent_ = nullptr;
  ContinuousTrajectory<Barycentric> const* trajectory_ = nullptr;
  ContinuousTrajectory<Barycentric>::Hint current_time_hint_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/celestial_body.hpp"
