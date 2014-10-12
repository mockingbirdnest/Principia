#pragma once

#include <memory>

#include "ksp_plugin/celestial.hpp"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"

using principia::physics::Body;
using principia::physics::Trajectory;
using principia::quantities::GravitationalParameter;

namespace principia {
namespace ksp_plugin {

// Represents a KSP |Vessel|.
template<typename Frame>
class Vessel {
 public:
  Vessel() = delete;
  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  ~Vessel() = default;

  // Constructs a vessel whose parent is initially |*parent|. |parent| must
  // not be null. No transfer of ownership.
  explicit Vessel(Celestial<Frame> const* parent);

  // A massless body.
  std::unique_ptr<Body<Frame> const> const body;
  // The parent body for the 2-body approximation. Not owning, must not be
  // null.
  Celestial<Frame> const* parent;
  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.
  std::unique_ptr<Trajectory<Frame>> history;
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
  Trajectory<Frame>* prolongation = nullptr;
  // Whether to keep the |Vessel| during the next call to |AdvanceTime|.
  bool keep = true;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/vessel_body.hpp"
