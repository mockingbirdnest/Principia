#pragma once

#include <string>

namespace principia {
namespace physics {

class SolarSystemFactory {
 public:
  void Initialize(std::string const& initial_state_filename,
                  std::string const& gravity_model_filename);

};

}  // namespace physics
}  // namespace principia

#include "physics/solar_system_factory_body.hpp"
