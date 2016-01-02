#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/dynamic_frame.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using quantities::Force;
using quantities::SpecificImpulse;
using physics::Frenet;
using geometry::Instant;

namespace ksp_plugin {

struct Burn {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
};

}  // namespace ksp_plugin
}  // namespace principia
