#pragma once

#include <string>

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

  int index(std::string const& name) const;

  static DegreesOfFreedom<Frame> MakeDegreesOfFreedom(
      serialization::InitialState::Body const& body);

  static std::unique_ptr<MassiveBody> MakeMassiveBody(
      serialization::GravityModel::Body const& body);

 private:
  std::vector<std::string> names_;
  std::unique_ptr<Ephemeris<Frame>> ephemeris_;
};

}  // namespace physics
}  // namespace principia

#include "physics/solar_system_factory_body.hpp"
