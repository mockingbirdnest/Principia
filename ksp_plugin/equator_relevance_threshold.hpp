#pragma once

#include "ksp_plugin/frames.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _equator_relevance_threshold {
namespace internal {

using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_quantities;

// Returns a distance from `body` that we consider is too far for the equator to
// be of interest.  Specifically, this distance is the maximum of
// - the semimajor axis of a supersynchronous orbit (1 orbit for 2 body
//   revolutions);
// - the distance at which a `Geopotential` with a tolerance of 0x1p-24 starts
//   damping dynamical oblateness;
Length EquatorRelevanceThreshold(RotatingBody<Barycentric> const& body);

}  // namespace internal

using internal::EquatorRelevanceThreshold;

}  // namespace _equator_relevance_threshold
}  // namespace ksp_plugin
}  // namespace principia
