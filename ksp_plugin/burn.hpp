
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
// mass.
struct Burn final {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::shared_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
  bool is_inertially_fixed;
};

}  // namespace internal_burn

using internal_burn::Burn;

}  // namespace ksp_plugin
}  // namespace principia
