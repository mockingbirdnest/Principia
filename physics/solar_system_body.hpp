#pragma once

#include "physics/solar_system.hpp"

#include <fstream>
#include <map>
#include <set>

#include "astronomy/frames.hpp"
#include "geometry/epoch.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/parser.hpp"
#include "quantities/si.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {

using geometry::Bivector;
using geometry::Instant;
using geometry::JulianDate;
using geometry::RadiusLatitudeLongitude;
using geometry::Vector;
using quantities::Length;
using quantities::ParseQuantity;
using quantities::Speed;
using quantities::si::Radian;
using quantities::si::Second;

namespace physics {

template<typename Frame>
void SolarSystem<Frame>::Initialize(std::string const& gravity_model_filename,
                                    std::string const& initial_state_filename) {
  // Parse the files.
  std::ifstream gravity_model_ifstream(gravity_model_filename);
  CHECK(gravity_model_ifstream.good());
  google::protobuf::io::IstreamInputStream gravity_model_zcs(
                                               &gravity_model_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&gravity_model_zcs,
                                            &gravity_model_));
  CHECK(gravity_model_.has_gravity_model());

  std::ifstream initial_state_ifstream(initial_state_filename);
  CHECK(initial_state_ifstream.good());
  google::protobuf::io::IstreamInputStream initial_state_zcs(
                                               &initial_state_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&initial_state_zcs,
                                            &initial_state_));
  CHECK(initial_state_.has_initial_state());

  // We don't support using a different frame that the one specified by the
  // instance.
  CHECK_EQ(Frame::tag, initial_state_.initial_state().frame());
  CHECK_EQ(Frame::tag, gravity_model_.gravity_model().frame());

  // Store the data in maps keyed by body name.
  for (auto const& body : gravity_model_.gravity_model().body()) {
    auto const inserted =
        gravity_model_map_.insert(std::make_pair(body.name(), &body));
    CHECK(inserted.second);
  }
  for (auto const& body : initial_state_.initial_state().body()) {
    auto const inserted =
        initial_state_map_.insert(std::make_pair(body.name(), &body));
    CHECK(inserted.second);
  }

  // Check that the maps are consistent.
  auto it1 = gravity_model_map_.begin();
  auto it2 = initial_state_map_.begin();
  for (; it1 != gravity_model_map_.end() && it2 != initial_state_map_.end();
       ++it1, ++it2) {
    CHECK_EQ(it1->first, it2->first);
    names_.push_back(it1->first);
  }
  CHECK(it1 == gravity_model_map_.end()) << it1->first;
  CHECK(it2 == initial_state_map_.end()) << it2->first;

  epoch_ = JulianDate(initial_state_.initial_state().epoch());

  // Call these two functions to parse all the data, so that errors are detected
  // at initialization.  Drop their results on the floor.
  MakeAllMassiveBodies();
  MakeAllDegreesOfFreedom();
}

template<typename Frame>
std::unique_ptr<Ephemeris<Frame>> SolarSystem<Frame>::MakeEphemeris(
    FixedStepSizeIntegrator<
        typename Ephemeris<Frame>::NewtonianMotionEquation> const&
        planetary_integrator,
    Time const& step,
    Length const& fitting_tolerance) {
  return std::make_unique<Ephemeris<Frame>>(MakeAllMassiveBodies(),
                                            MakeAllDegreesOfFreedom(), 
                                            epoch_,
                                            planetary_integrator,
                                            step,
                                            fitting_tolerance);
}

template<typename Frame>
Instant const& SolarSystem<Frame>::epoch() const {
  return epoch_;
}

template<typename Frame>
std::vector<std::string> const& SolarSystem<Frame>::names() const {
  return names_;
}

template<typename Frame>
int SolarSystem<Frame>::index(std::string const& name) const {
  auto const it = std::equal_range(names_.begin(), names_.end(), name);
  return it.first - names_.begin();
}

template<typename Frame>
MassiveBody const& SolarSystem<Frame>::massive_body(
    Ephemeris<Frame> const & ephemeris,
    std::string const & name) const {
  return *ephemeris.bodies()[index(name)];
}

template<typename Frame>
ContinuousTrajectory<Frame> const& SolarSystem<Frame>::trajectory(
    Ephemeris<Frame> const & ephemeris,
    std::string const & name) const {
  MassiveBody const* const body = ephemeris.bodies()[index(name)];
  return *ephemeris.trajectory(body);
}

template<typename Frame>
serialization::InitialState::Body const&
SolarSystem<Frame>::initial_state(std::string const& name) const {
  return *FindOrDie(initial_state_map_, name);
}

template<typename Frame>
serialization::GravityModel::Body const&
SolarSystem<Frame>::gravity_model(std::string const& name) const {
  return *FindOrDie(gravity_model_map_, name);
}

template<typename Frame>
DegreesOfFreedom<Frame> SolarSystem<Frame>::MakeDegreesOfFreedom(
    serialization::InitialState::Body const& body) {
  Position<Frame> const
      position = Frame::origin + 
                 Vector<Length, Frame>({ParseQuantity<Length>(body.x()),
                                        ParseQuantity<Length>(body.y()),
                                        ParseQuantity<Length>(body.z())});
  Velocity<Frame> const
      velocity(Vector<Speed, Frame>({ParseQuantity<Speed>(body.vx()),
                                     ParseQuantity<Speed>(body.vy()),
                                     ParseQuantity<Speed>(body.vz())}));
  return DegreesOfFreedom<Frame>(position, velocity);
}

template<typename Frame>
std::unique_ptr<MassiveBody> SolarSystem<Frame>::MakeMassiveBody(
    serialization::GravityModel::Body const& body) {
  CHECK(body.has_gravitational_parameter());
  CHECK_EQ(body.has_j2(), body.has_reference_radius());
  CHECK_EQ(body.has_axis_declination(), body.has_axis_right_ascension());
  MassiveBody::Parameters massive_body_parameters(
                              ParseQuantity<GravitationalParameter>(
                                  body.gravitational_parameter()));
  if (body.has_axis_declination()) {
    // TODO(phl): Parse the additional parameters.
    RotatingBody<Frame>::Parameters
        rotating_body_parameters(
            0 * Radian,
            Instant(),
            Bivector<double, Frame>(
                RadiusLatitudeLongitude(
                    1.0,
                    ParseQuantity<Angle>(body.axis_declination()),
                    ParseQuantity<Angle>(body.axis_right_ascension())).
                ToCartesian()) * Radian / Second);
    if (body.has_j2()) {
      OblateBody<Frame>::Parameters oblate_body_parameters(
          ParseQuantity<double>(body.j2()),
          ParseQuantity<Length>(body.reference_radius()));
      return std::make_unique<OblateBody<Frame>>(massive_body_parameters,
                                                 rotating_body_parameters,
                                                 oblate_body_parameters);
    } else {
      return std::make_unique<RotatingBody<Frame>>(massive_body_parameters,
                                                   rotating_body_parameters);
    }
  } else {
    return std::make_unique<MassiveBody>(massive_body_parameters);
  }
}

template<typename Frame>
std::vector<not_null<std::unique_ptr<MassiveBody const>>>
SolarSystem<Frame>::MakeAllMassiveBodies() {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  for (auto const& pair : gravity_model_map_) {
    std::string const& name = pair.first;
    serialization::GravityModel::Body const* const body = pair.second;
    CHECK(body->has_gravitational_parameter());
    CHECK_EQ(body->has_j2(), body->has_reference_radius());
    CHECK_EQ(body->has_axis_declination(), body->has_axis_right_ascension());
    bodies.emplace_back(MakeMassiveBody(*body));
  }
  return bodies;
}

template<typename Frame>
std::vector<DegreesOfFreedom<Frame>>
SolarSystem<Frame>::MakeAllDegreesOfFreedom() {
  std::vector<DegreesOfFreedom<Frame>> degrees_of_freedom;
  for (auto const& pair : initial_state_map_) {
    std::string const& name = pair.first;
    serialization::InitialState::Body const* const body = pair.second;
    degrees_of_freedom.push_back(MakeDegreesOfFreedom(*body));
  }
  return degrees_of_freedom;
}

}  // namespace physics
}  // namespace principia
