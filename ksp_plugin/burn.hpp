#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "physics/dynamic_frame.hpp"

namespace principia {
namespace ksp_plugin {

struct Burn {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> ?v;
};

}  // namespace ksp_plugin
}  // namespace principia
