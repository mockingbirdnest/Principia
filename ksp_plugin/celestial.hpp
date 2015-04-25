#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/mobile_interface.hpp"
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::not_null;
using physics::Body;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using quantities::GravitationalParameter;

namespace ksp_plugin {

// Represents a KSP |CelestialBody|.
class Celestial : public MobileInterface {
 public:
  explicit Celestial(not_null<std::unique_ptr<MassiveBody const>> body);
  Celestial(Celestial const&) = delete;
  Celestial(Celestial&&) = delete;
  ~Celestial() = default;

  // True if, and only if, |history_| is not null.
  bool is_initialized() const;

  MassiveBody const& body() const;
  bool has_parent() const;
  Celestial const* parent() const;  // Null for the Sun.
  void set_parent(not_null<Celestial const*> const parent);

  // Both accessors require |is_initialized()|.
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

  // Creates a |history_| for this body and appends a point with the given
  // |time| and |degrees_of_freedom|.  Then forks a |prolongation_| at |time|.
  // The celestial |is_initialized()| after the call.
  void CreateHistoryAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Deletes the |prolongation_| and forks a new one at |time|.
  void ResetProlongation(Instant const& time);

  // Creates a |prediction_| forked at the end of the |prolongation_|.
  // Requires |is_initialized()| and |!has_prediction()|.
  void ForkPrediction();

  // Deletes the |prediction_|.
  void DeletePrediction();

  // The celestial must satisfy |is_initialized()|.
  void WriteToMessage(not_null<serialization::Celestial*> const message) const;
  // NOTE(egg): This should return a |not_null|, but we can't do that until
  // |not_null<std::unique_ptr<T>>| is convertible to |std::unique_ptr<T>|, and
  // that requires a VS 2015 feature (rvalue references for |*this|).
  static std::unique_ptr<Celestial> ReadFromMessage(
      serialization::Celestial const& message);

 private:
  not_null<std::unique_ptr<MassiveBody const>> const body_;
  // The parent body for the 2-body approximation. Not owning, must only
  // be null for the sun.
  Celestial const* parent_ = nullptr;
  // The past and present trajectory of the body. It ends at |HistoryTime()|.
  std::unique_ptr<Trajectory<Barycentric>> history_;
  // A child trajectory of |*history|. It is forked at |history->last_time()|
  // and continues it until |current_time_|. It is computed with a
  // non-constant timestep, which breaks symplecticity. |history| is advanced
  // with a constant timestep as soon as possible, and |prolongation| is then
  // restarted from this new end of |history|.
  // Not owning.
  Trajectory<Barycentric>* prolongation_ = nullptr;
  // A child trajectory of |*prolongation_|.
  Trajectory<Barycentric>* prediction_ = nullptr;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/celestial_body.hpp"
