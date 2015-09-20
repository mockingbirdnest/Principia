#pragma once

#include "physics/solar_system_factory.hpp"

#include <fstream>
#include <map>
#include <set>

#include "astronomy/frames.hpp"
#include "geometry/grassmann.hpp"
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
using geometry::RadiusLatitudeLongitude;
using geometry::Vector;
using quantities::Length;
using quantities::ParseQuantity;
using quantities::Speed;
using quantities::si::Radian;
using quantities::si::Second;

namespace physics {

namespace {

template<typename Frame>
DegreesOfFreedom<Frame> MakeDegreesOfFreedom(
    serialization::InitialState::Body const& body);

template<typename Frame>
std::unique_ptr<MassiveBody> MakeMassiveBody(
    serialization::GravityModel::Body const& body);

template<typename Frame>
DegreesOfFreedom<Frame> MakeDegreesOfFreedom(
    serialization::InitialState::Body const& body) {
  Position<Frame> const
      position(Vector<Length, Frame>({ParseQuantity<Length>(body.x()),
                                      ParseQuantity<Length>(body.y()),
                                      ParseQuantity<Length>(body.z())}));
  Velocity<Frame> const
      velocity(Vector<Speed, Frame>({ParseQuantity<Speed>(body.vx()),
                                     ParseQuantity<Speed>(body.vy()),
                                     ParseQuantity<Speed>(body.vz())}));
  return DegreesOfFreedom<Frame>(position, velocity);
}

template<typename Frame>
std::unique_ptr<MassiveBody> MakeMassiveBody(
    serialization::GravityModel::Body const& body) {
  CHECK(body.has_gravitational_parameter());
  CHECK_EQ(body.has_j2(), body.has_reference_radius());
  CHECK_EQ(body.has_axis_declination(), body.has_axis_right_ascension());
  MassiveBody::Parameters massive_body_parameters(
                              ParseQuantity<GravitationalParameter>(
                                  body.gravitational_parameter()));
  if (body.has_axis_declination()) {
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

}  // namespace

void SolarSystemFactory::Initialize(std::string const& gravity_model_filename,
                                    std::string const& initial_state_filename) {
  // Parse the files.
  serialization::SolarSystemFile gravity_model;
  std::ifstream gravity_model_ifstream(gravity_model_filename);
  CHECK(gravity_model_ifstream.good());
  google::protobuf::io::IstreamInputStream gravity_model_zcs(
                                               &gravity_model_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&gravity_model_zcs,
                                            &gravity_model));
  CHECK(gravity_model.has_gravity_model());

  serialization::SolarSystemFile initial_state;
  std::ifstream initial_state_ifstream(initial_state_filename);
  CHECK(initial_state_ifstream.good());
  google::protobuf::io::IstreamInputStream initial_state_zcs(
                                               &initial_state_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&initial_state_zcs,
                                            &initial_state));
  CHECK(initial_state.has_initial_state());

  // We don't support using different frames in different files.
  CHECK_EQ(gravity_model.gravity_model().frame(),
           initial_state.initial_state().frame());
  serialization::Frame::SolarSystemTag const frame =
      gravity_model.initial_state().frame();

  // Store the data in maps keyed by body name.
  std::map<std::string,
           serialization::GravityModel::Body const*> gravity_model_map;
  for (auto const& body : gravity_model.gravity_model().body()) {
    auto const inserted =
        gravity_model_map.insert(std::make_pair(body.name(), &body));
    CHECK(inserted.second);
  }
  std::map<std::string,
           serialization::InitialState::Body const*> initial_state_map;
  for (auto const& body : initial_state.initial_state().body()) {
    auto const inserted =
        initial_state_map.insert(std::make_pair(body.name(), &body));
    CHECK(inserted.second);
  }

  // Check that the maps are consistent.
  auto it1 = gravity_model_map.begin();
  auto it2 = initial_state_map.begin();
  for (; it1 != gravity_model_map.end() && it2 != initial_state_map.end();
       ++it1, ++it2) {
    CHECK_EQ(it1->first, it2->first);
  }
  CHECK(it1 == gravity_model_map.end()) << it1->first;
  CHECK(it2 == initial_state_map.end()) << it2->first;

  // Build bodies.
  for (auto const& pair : gravity_model_map) {
    std::string const& name = pair.first;
    serialization::GravityModel::Body const* const body = pair.second;
    std::unique_ptr<MassiveBody> massive_body;
    CHECK(body->has_gravitational_parameter());
    CHECK_EQ(body->has_j2(), body->has_reference_radius());
    CHECK_EQ(body->has_axis_declination(), body->has_axis_right_ascension());
    switch (frame) {
      case serialization::Frame::ICRF_J2000_ECLIPTIC:
        massive_body = MakeMassiveBody<astronomy::ICRFJ2000Ecliptic>(*body);
        break;
      case serialization::Frame::ICRF_J2000_EQUATOR:
        massive_body = MakeMassiveBody<astronomy::ICRFJ2000Equator>(*body);
        break;
    }
  }

  // Build degrees of freedom.
  for (auto const& pair : initial_state_map) {
    std::string const& name = pair.first;
    serialization::InitialState::Body const* const body = pair.second;
    switch (frame) {
      case serialization::Frame::ICRF_J2000_ECLIPTIC: {
        auto const degrees_of_freedom =
            MakeDegreesOfFreedom<astronomy::ICRFJ2000Ecliptic>(*body);
        break;
      }
      case serialization::Frame::ICRF_J2000_EQUATOR: {
        auto const degrees_of_freedom =
            MakeDegreesOfFreedom<astronomy::ICRFJ2000Equator>(*body);
        break;
      }
    }
  }
}

}  // namespace physics
}  // namespace principia
