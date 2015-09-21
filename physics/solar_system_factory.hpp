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
class SolarSystemFactory {
 public:
  void Initialize(std::string const& gravity_model_filename,
                  std::string const& initial_state_filename);

  std::unique_ptr<Ephemeris<Frame>> MakeEphemeris(
      FixedStepSizeIntegrator<
          typename Ephemeris<Frame>::NewtonianMotionEquation> const&
          planetary_integrator,
      Time const& step,
      Length const& fitting_tolerance);

  int index(std::string const& name) const;
  MassiveBody const& massive_body(Ephemeris<Frame> const& ephemeris,
                                  std::string const& name) const;
  ContinuousTrajectory<Frame> const& trajectory(
      Ephemeris<Frame> const& ephemeris,
      std::string const& name) const;

  static DegreesOfFreedom<Frame> MakeDegreesOfFreedom(
      serialization::InitialState::Body const& body);

  static std::unique_ptr<MassiveBody> MakeMassiveBody(
      serialization::GravityModel::Body const& body);

 private:
  Instant epoch_;
  std::vector<std::string> names_;
  std::map<std::string,
           serialization::GravityModel::Body const*> gravity_model_map_;
  std::map<std::string,
           serialization::InitialState::Body const*> initial_state_map_;
};

}  // namespace physics
}  // namespace principia

#include "physics/solar_system_factory_body.hpp"
