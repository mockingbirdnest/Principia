#pragma once

#include <memory>
#include <utility>

#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"

using principia::physics::Body;
using principia::physics::DegreesOfFreedom;
using principia::physics::Trajectory;
using principia::quantities::GravitationalParameter;

namespace principia {
namespace ksp_plugin {

// Represents a KSP |CelestialBody|.
template<typename Frame>
class Celestial {
 public:
  Celestial(Celestial const&) = delete;
  Celestial(Celestial&&) = delete;
  ~Celestial() = default;

  template<typename... Args>
  explicit Celestial(Args&&... args);  // NOLINT(build/c++11)

  Body<Frame> const& body() const;
  bool has_parent() const;
  Celestial const& parent() const;
  Trajectory<Frame> const& history() const;
  Trajectory<Frame> const& prolongation() const;

  Trajectory<Frame>* mutable_history();
  Trajectory<Frame>* mutable_prolongation();
  void set_parent(Celestial const* parent);

  // Creates a |history_| for this body and appends a point with the given
  // |time| and |degrees_of_freedom|.  Then forks a |prolongation_| at |time|.
  void AppendAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Deletes the |prolongation_| and forks a new one at |time|.
  void ResetProlongation(Instant const& time);

 private:
  std::unique_ptr<Body<Frame> const> const body_;
  // The parent body for the 2-body approximation. Not owning, must only
  // be null for the sun.
  Celestial const* parent_ = nullptr;
  // The past and present trajectory of the body. It ends at |HistoryTime()|.
  std::unique_ptr<Trajectory<Frame>> history_;
  // A child trajectory of |*history|. It is forked at |history->last_time()|
  // and continues it until |current_time_|. It is computed with a
  // non-constant timestep, which breaks symplecticity. |history| is advanced
  // with a constant timestep as soon as possible, and |prolongation| is then
  // restarted from this new end of |history|.
  // Not owning, not null.
  Trajectory<Frame>* prolongation_ = nullptr;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/celestial_body.hpp"
