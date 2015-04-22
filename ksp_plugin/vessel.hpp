#pragma once

#include <memory>

#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/mobile_interface.hpp"
#include "ksp_plugin/vessel.hpp"
#include "ksp_plugin/part.hpp"
#include "physics/massless_body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using physics::MasslessBody;
using physics::Trajectory;
using quantities::GravitationalParameter;

namespace ksp_plugin {

// Represents a KSP |Vessel|.
class Vessel : public MobileInterface {
 public:
  Vessel() = delete;
  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  Vessel& operator=(Vessel const&) = delete;
  Vessel& operator=(Vessel&&) = delete;
  ~Vessel() = default;

  // Constructs a vessel whose parent is initially |*parent|.  No transfer of
  // ownership.
  explicit Vessel(not_null<Celestial const*> const parent);

  // True if, and only if, |history_| is not null.
  bool is_synchronized() const;
  // True if, and only if, |prolongation_| is not null, i.e., if either
  // |CreateProlongation| or |CreateHistoryAndForkProlongation| was called at
  // some point.
  bool is_initialized() const;

  not_null<Celestial const*> parent() const;
  void set_parent(not_null<Celestial const*> const parent);

  // Both accessors require |is_synchronized()|.
  Trajectory<Barycentric> const& history() const override;
  not_null<Trajectory<Barycentric>*> mutable_history() override;

  // Both accessors require |is_initialized()|.
  Trajectory<Barycentric> const& prolongation() const override;
  not_null<Trajectory<Barycentric>*> mutable_prolongation() override;

  // Both accessors require |is_initialized()|.  In addition the first one
  // requires |has_prediction()|.
  Trajectory<Barycentric> const& prediction() const override;
  Trajectory<Barycentric>* mutable_prediction() override;
  bool has_prediction() const override;

  // Creates an |owned_prolongation_| for this vessel and appends a point with
  // the given |time| and |degrees_of_freedom|.  The vessel must not satisfy
  // |is_initialized()| nor |is_synchronized()|, |owned_prolongation_| must be
  // null.  The vessel |is_initialized()|, but does not satisfy
  // |is_synchronized()|, after the call.
  void CreateProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Creates a |history_| for this vessel and appends a point with the
  // given |time| and |degrees_of_freedom|, then forks a |prolongation_| at
  // |time|.  Nulls |owned_prolongation_|.  The vessel must not satisfy
  // |is_synchronized()|.  |*owned_prolongation_| is destroyed *after*
  // |history_| has been constructed.
  // The vessel |is_synchronized()| and |is_initialized()| after the call.
  void CreateHistoryAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Deletes the |prolongation_| and forks a new one at |time|.
  // The vessel must satisfy |is_synchronized()| and |is_initialized()|,
  // |owned_prolongation_| must be null.
  void ResetProlongation(Instant const& time);

  // Creates a |prediction_| forked at the end of the |prolongation_|.
  // Requires |is_initialized()| and |!has_prediction()|.
  void ForkPrediction();

  // Deletes the |prediction_|.
  void DeletePrediction();

  // The vessel must satisfy |is_initialized()|.
  void WriteToMessage(not_null<serialization::Vessel*> const message) const;
  // NOTE(egg): This should return a |not_null|, but we can't do that until
  // |not_null<std::unique_ptr<T>>| is convertible to |std::unique_ptr<T>|, and
  // that requires a VS 2015 feature (rvalue references for |*this|).
  static std::unique_ptr<Vessel> ReadFromMessage(
      serialization::Vessel const& message,
      not_null<Celestial const*> const parent);

 private:
  MasslessBody const body_;
  // The parent body for the 2-body approximation. Not owning.
  not_null<Celestial const*> parent_;
  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.  It is advanced with a constant time step.
  std::unique_ptr<Trajectory<Barycentric>> history_;
  // Most of the time, this is a child trajectory of |*history_|. It is forked
  // at |history_->last_time()| and continues until |current_time_|. It is
  // computed with a non-constant timestep, which breaks symplecticity.
  // If |history_| is null, this points to |owned_prolongation_| instead.
  // Not owning.
  Trajectory<Barycentric>* prolongation_ = nullptr;
  // When the vessel is added, before it is synchonized with the other vessels
  // and celestials, there is no |history_|.  The prolongation is directly owned
  // during that time.  Null if, and only if, |history_| is not null.
  std::unique_ptr<Trajectory<Barycentric>> owned_prolongation_;
  // A child trajectory of |*prolongation_|.
  Trajectory<Barycentric>* prediction_ = nullptr;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/vessel_body.hpp"
