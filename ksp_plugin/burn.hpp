
#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_burn {

using base::not_null;
using geometry::Instant;
using geometry::Velocity;
using physics::Frenet;
using quantities::Force;
using quantities::Mass;
using quantities::SpecificImpulse;

// Parameters for constructing a |NavigationManœuvre|, excluding the initial
// mass.  This owns a |NavigationFrame| and is therefore not copyable.
struct Burn final {
  // TODO(egg): we won't need that constructor once the compiler supports C++14.
  Burn(Force const thrust,
       SpecificImpulse const specific_impulse,
       not_null<std::unique_ptr<NavigationFrame const>> frame,
       Instant const initial_time,
       Velocity<Frenet<Navigation>> const& Δv,
       bool const is_inertially_fixed)
      : thrust(thrust),
        specific_impulse(specific_impulse),
        frame(std::move(frame)),
        Δv(Δv),
        is_inertially_fixed(is_inertially_fixed) {}

  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
  bool is_inertially_fixed{};
};

NavigationManœuvre MakeNavigationManœuvre(Burn burn, Mass const& initial_mass);

}  // namespace internal_burn

using internal_burn::Burn;
using internal_burn::MakeNavigationManœuvre;

}  // namespace ksp_plugin
}  // namespace principia
