#pragma once

#include <memory>

#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"

using principia::physics::Body;
using principia::physics::Trajectory;
using principia::quantities::GravitationalParameter;

namespace principia {
namespace ksp_plugin {

// Represents a KSP |CelestialBody|.
template<typename Frame>
struct Celestial {
  Celestial() = delete;
  Celestial(Celestial const&) = delete;
  Celestial(Celestial&&) = delete;
  ~Celestial() = default;

  explicit Celestial(GravitationalParameter const& gravitational_parameter)
    : body(new Body<Frame>(gravitational_parameter)) {}
  std::unique_ptr<Body<Frame> const> const body;
  // The parent body for the 2-body approximation. Not owning, must only
  // be null for the sun.
  Celestial const* parent = nullptr;
  // The past and present trajectory of the body. It ends at |HistoryTime()|.
  std::unique_ptr<Trajectory<Frame>> history;
  // A child trajectory of |*history|. It is forked at |history->last_time()|
  // and continues it until |current_time_|. It is computed with a
  // non-constant timestep, which breaks symplecticity. |history| is advanced
  // with a constant timestep as soon as possible, and |prolongation| is then
  // restarted from this new end of |history|.
  // Not owning, not null.
  Trajectory<Frame>* prolongation = nullptr;
};

}  // namespace ksp_plugin
}  // namespace principia
