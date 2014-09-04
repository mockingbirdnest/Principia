#pragma once

#include <map>
#include <string>

#include "physics/body.hpp"

namespace principia {
namespace ksp_plugin {

class Plugin {
 public:
 private:
  std::map<std::string, physics::Body> vessels_;
  std::map<int, physics::Body> celestial_bodies_;
};

}  // namespace ksp_plugin
}  // namespace principia
