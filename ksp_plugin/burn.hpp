#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/dynamic_frame.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Instant;
using geometry::Velocity;
using physics::DynamicFrame;
using physics::Frenet;
using quantities::Force;
using quantities::SpecificImpulse;

namespace ksp_plugin {

//TODO(phl): Not nice that this is here.
using NavigationFrame = DynamicFrame<Barycentric, Navigation>;
using NavigationManœuvre = Manœuvre<Barycentric, Navigation>;

// Parameters for constructing a |NavigationManœuvre|, excluding the initial
// mass.  This owns a |NavigationFrame| and is therefore not copyable.
struct Burn {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
};

NavigationManœuvre MakeNavigationManœuvre(Burn burn, Mass const& initial_mass);

}  // namespace ksp_plugin
}  // namespace principia
