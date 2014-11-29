#pragma once

#include <memory>

#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/vessel.hpp"
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

  Trajectory<Barycentric>* mutable_history();
  Trajectory<Barycentric>* mutable_prolongation();
  void set_parent(Celestial const* parent);

  // Creates an |owned_prolongation_| for this vessel and appends a point with
  // the given |time| and |degrees_of_freedom|.  The vessel must have neither a
  // |history_|, an |owned_prolongation_| nor a |prolongation_|.
  void CreateProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Creates a |history_| for this vessel and appends a point with  with the
  // given |time| and |degrees_of_freedom|, then forks a |prolongation_| at
  // |time|.  Nulls |owned_prolongation_|.  The vessel must not have a
  // |history_|. |*owned_prolongation_| is destroyed *after* |history_| has been
  // constructed.
  void CreateHistoryAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Deletes the |prolongation_| and forks a new one at |time|.
  // |history_| and |prolongation_| must not be null, |owned_prolongation_|
  // must be null.
  void ResetProlongation(Instant const& time);

 private:
  std::unique_ptr<MasslessBody const> const body_;
  // The parent body for the 2-body approximation. Not owning, must not be
  // null.
  Celestial const* parent_;
  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.  It is advanced with a constant time step.
  std::unique_ptr<Trajectory<Barycentric>> history_;
  // Most of the time, this is a child trajectory of |*history_|. It is forked
  // at/ |history_->last_time()| and continues it until |current_time_|. It is
  // computed with a non-constant timestep, which breaks symplecticity.
  // If |history_| is null, this points to |owned_prolongation_| instead.
  // Not owning.
  Trajectory<Barycentric>* prolongation_ = nullptr;
  // When the vessel is added, before it is synchonized with the other vessels
  // and celestials, there is no |history|.  The prolongation is directly owned
  // during that time.  Null if, and only if, |history_| is not null.
  std::unique_ptr<Trajectory<Barycentric>> owned_prolongation_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/vessel_body.hpp"
