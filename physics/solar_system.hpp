
#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "base/not_null.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/body.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/hierarchical_system.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {
namespace physics {
namespace internal_solar_system {

using base::not_null;
using geometry::Instant;
using quantities::GravitationalParameter;
using quantities::Length;

template<typename Frame>
class SolarSystem final {
 public:
  // Constructs a solar system from the given files, which must contain text
  // format for SolarSystemFile protocol buffers.
  SolarSystem(std::filesystem::path const& gravity_model_filename,
              std::filesystem::path const& initial_state_filename,
              bool ignore_frame = false);

  // Construct a solar system from the given messages.
  SolarSystem(serialization::GravityModel const& gravity_model,
              serialization::InitialState const& initial_state,
              bool ignore_frame = false);

  SolarSystem(SolarSystem const& other);
  SolarSystem& operator=(const SolarSystem& other);

  // Constructs an ephemeris for this object using the specified parameters.
  // The bodies and initial state are constructed from the data passed to
  // |Initialize|.
  not_null<std::unique_ptr<Ephemeris<Frame>>> MakeEphemeris(
      Length const& fitting_tolerance,
      typename Ephemeris<Frame>::FixedStepParameters const& parameters) const;

  std::vector<not_null<std::unique_ptr<MassiveBody const>>>
  MakeAllMassiveBodies() const;

  // The time origin for the initial state.
  Instant const& epoch() const;
  std::string const& epoch_literal() const;

  // The names of the bodies, sorted alphabetically.
  std::vector<std::string> const& names() const;

  // The index of the body named |name| in the vector |names()| and in the
  // bodies of the ephemeris.
  int index(std::string const& name) const;

  // The initial state of the body named |name|.
  DegreesOfFreedom<Frame> degrees_of_freedom(std::string const& name) const;

  // The gravitational parameter of the body named |name|.
  GravitationalParameter gravitational_parameter(std::string const& name) const;

  // The mean radius of the body named |name|.
  Length mean_radius(std::string const& name) const;

  // The |MassiveBody| for the body named |name|, extracted from the given
  // |ephemeris|.
  not_null<MassiveBody const*> massive_body(Ephemeris<Frame> const& ephemeris,
                                            std::string const& name) const;

  // Same as above, but checks that the body is a rotating body.
  not_null<RotatingBody<Frame> const*> rotating_body(
      Ephemeris<Frame> const& ephemeris,
      std::string const& name) const;

  // The |ContinuousTrajectory| for the body named |name|, extracted from the
  // given |ephemeris|.
  ContinuousTrajectory<Frame> const& trajectory(
      Ephemeris<Frame> const& ephemeris,
      std::string const& name) const;

  // The configuration protocol buffers for the body named |name|.
  serialization::GravityModel::Body const& gravity_model_message(
      std::string const& name) const;
  bool has_cartesian_initial_state_message(std::string const& name) const;
  serialization::InitialState::Cartesian::Body const&
  cartesian_initial_state_message(std::string const& name) const;
  bool has_keplerian_initial_state_message(std::string const& name) const;
  serialization::InitialState::Keplerian::Body const&
  keplerian_initial_state_message(std::string const& name) const;

  // Factory functions for converting configuration protocol buffers into
  // structured objects.
  static DegreesOfFreedom<Frame> MakeDegreesOfFreedom(
      serialization::InitialState::Cartesian::Body const& body);
  static KeplerianElements<Frame> MakeKeplerianElements(
      serialization::InitialState::Keplerian::Body::Elements const& elements);
  static not_null<std::unique_ptr<MassiveBody>> MakeMassiveBody(
      serialization::GravityModel::Body const& body);
  static not_null<std::unique_ptr<RotatingBody<Frame>>> MakeRotatingBody(
      serialization::GravityModel::Body const& body);
  static not_null<std::unique_ptr<OblateBody<Frame>>> MakeOblateBody(
      serialization::GravityModel::Body const& body);

  // Constructs a hierarchical system for a solar system defined by keplerian
  // elements.
  not_null<std::unique_ptr<HierarchicalSystem<Frame>>> MakeHierarchicalSystem()
      const;

  // Utilities for patching the internal protocol buffers after initialization.
  void LimitOblatenessToDegree(std::string const& name, int max_degree);
  void LimitOblatenessToZonal(std::string const& name);
  void RemoveMassiveBody(std::string const& name);
  void ReplaceElements(std::string const& name,
                       KeplerianElements<Frame> const& elements);

 private:
  // Fails if the given |body| doesn't have a consistent set of fields.
  static void Check(serialization::GravityModel::Body const& body);

  // Factory functions for the parameter classes of the bodies.
  static not_null<std::unique_ptr<MassiveBody::Parameters>>
  MakeMassiveBodyParameters(serialization::GravityModel::Body const& body);
  static not_null<std::unique_ptr<typename RotatingBody<Frame>::Parameters>>
  MakeRotatingBodyParameters(serialization::GravityModel::Body const& body);
  static not_null<std::unique_ptr<typename OblateBody<Frame>::Parameters>>
  MakeOblateBodyParameters(serialization::GravityModel::Body const& body);

  std::vector<DegreesOfFreedom<Frame>> MakeAllDegreesOfFreedom() const;

  // If a frame is specified in a message it must match the frame of this
  // instance.  Otherwise the frame of the instance is used.  This is convenient
  // for tests.
  template<typename Message>
  static void CheckFrame(Message const& message);

  // Note that the maps below hold pointers into these protocol buffers.
  serialization::GravityModel gravity_model_;
  serialization::InitialState initial_state_;

  Instant epoch_;
  std::vector<std::string> names_;
  std::map<std::string,
           serialization::GravityModel::Body*> gravity_model_map_;
  std::map<std::string,
           serialization::InitialState::Cartesian::Body const*>
      cartesian_initial_state_map_;
  std::map<std::string,
           serialization::InitialState::Keplerian::Body*>
      keplerian_initial_state_map_;
};

}  // namespace internal_solar_system

using internal_solar_system::SolarSystem;

}  // namespace physics
}  // namespace principia

#include "physics/solar_system_body.hpp"
