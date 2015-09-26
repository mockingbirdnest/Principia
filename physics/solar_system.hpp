#pragma once

#include <map>
#include <string>
#include <vector>

#include "integrators/ordinary_differential_equations.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {
namespace physics {

template<typename Frame>
class SolarSystem {
 public:
  // Initializes this object from the given files, which must contain text
  // format for SolarSystemFile protocol buffers.
  void Initialize(std::string const& gravity_model_filename,
                  std::string const& initial_state_filename);

  // Constructs an ephemeris for this object using the specified parameters.
  // The bodies and initial state are constructed from the data passed to
  // |Initialize|.
  std::unique_ptr<Ephemeris<Frame>> MakeEphemeris(
      FixedStepSizeIntegrator<
          typename Ephemeris<Frame>::NewtonianMotionEquation> const&
          planetary_integrator,
      Time const& step,
      Length const& fitting_tolerance);

  // The time origin for the initial state.
  Instant const& epoch() const;

  // The names of the bodies, sorted alphabetically.
  std::vector<std::string> const& names() const;

  // The index of the body named |name| in the vector |names()| and in the
  // bodies of the ephemeris.
  int index(std::string const& name) const;

  // The |MassiveBody| for the body named |name|, extracted from the given
  // |ephemeris|.
  MassiveBody const& massive_body(Ephemeris<Frame> const& ephemeris,
                                  std::string const& name) const;

  // The |ContinuousTrajectory| for the body named |name|, extracted from the
  // given |ephemeris|.
  ContinuousTrajectory<Frame> const& trajectory(
      Ephemeris<Frame> const& ephemeris,
      std::string const& name) const;

  // The configuration protocol buffers for the body named |name|.
  serialization::InitialState::Body const& initial_state(
      std::string const& name) const;
  serialization::GravityModel::Body const& gravity_model(
      std::string const& name) const;

  // Factory functions for converting configuration protocol buffers into
  // structured objects.
  static DegreesOfFreedom<Frame> MakeDegreesOfFreedom(
      serialization::InitialState::Body const& body);
  static std::unique_ptr<MassiveBody> MakeMassiveBody(
      serialization::GravityModel::Body const& body);

  // Utilities for patching the internal protocol buffers after initialization.
  // Should only be used in tests.
  void RemoveMassiveBody(std::string const& name);
  void RemoveOblateness(std::string const& name);

 private:
  std::vector<not_null<std::unique_ptr<MassiveBody const>>>
      MakeAllMassiveBodies();
  std::vector<DegreesOfFreedom<Frame>> MakeAllDegreesOfFreedom();

  // Note that the maps below hold pointers into these protocol buffers.
  serialization::SolarSystemFile gravity_model_;
  serialization::SolarSystemFile initial_state_;

  Instant epoch_;
  std::vector<std::string> names_;
  std::map<std::string,
           serialization::InitialState::Body const*> initial_state_map_;
  std::map<std::string,
           serialization::GravityModel::Body*> gravity_model_map_;
};

}  // namespace physics
}  // namespace principia

#include "physics/solar_system_body.hpp"
