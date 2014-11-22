#pragma once

#include <memory>

#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/vessel.hpp"
#include "ksp_plugin/part.hpp"
#include "physics/massless_body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"

using principia::physics::MasslessBody;
using principia::physics::Trajectory;
using principia::quantities::GravitationalParameter;

namespace principia {
namespace ksp_plugin {

// Represents a KSP |Vessel|.
class Vessel {
 public:
  Vessel() = delete;
  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  ~Vessel() = default;

  // Constructs a vessel whose parent is initially |*parent|. |parent| must
  // not be null. No transfer of ownership.
  explicit Vessel(Celestial const* parent);

  bool has_history() const;
  bool has_prolongation() const;
  Celestial const& parent() const;
  Trajectory<Barycentric> const& history() const;
  Trajectory<Barycentric> const& prolongation() const;
  Trajectory<Barycentric> const& prolongation_or_history() const;

  Trajectory<Barycentric>* mutable_history();
  Trajectory<Barycentric>* mutable_prolongation();
  void set_parent(Celestial const* parent);

  // Creates a |history_| for this body and appends a point with the given
  // |time| and |degrees_of_freedom|.
  void CreateHistory(Instant const& time,
                     DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Deletes the |prolongation_| if there is one and forks a new one at |time|.
  void ResetProlongation(Instant const& time);

 private:
  std::unique_ptr<MasslessBody const> const body_;
  // The parent body for the 2-body approximation. Not owning, must not be
  // null.
  Celestial const* parent_;
  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.
  std::unique_ptr<Trajectory<Barycentric>> history_;
  // A child trajectory of |*history|. It is forked at |history->last_time()|
  // and continues it until |current_time_|. It is computed with a
  // non-constant timestep, which breaks symplecticity. |history| is advanced
  // with a constant timestep as soon as possible, and |prolongation| is then
  // restarted from this new end of |history|.
  // Not owning, is null when the vessel is added and becomes non-null when
  // |history| is next advanced for all vessels and celestials. In the
  // meantime, |history| is advanced with small, non-constant timesteps to
  // catch up with the synchronous constant-timestep integration.
  // |this| is in |new_vessels_| if and only if |prolongation| is null.
  Trajectory<Barycentric>* prolongation_ = nullptr;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/vessel_body.hpp"
