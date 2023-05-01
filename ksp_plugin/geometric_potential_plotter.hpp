#pragma once

#include "base/jthread.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/equipotential.hpp"

namespace principia {
namespace ksp_plugin {
namespace _geometric_potential_plotter {
namespace internal {

using namespace principia::base::_jthread;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_equipotential;
using namespace principia::physics::_massive_body;
using namespace principia::ksp_plugin::_frames;

class GeometricPotentialPlotter {
 public:
  struct Parameters {
    not_null<MassiveBody const*> primary;
    not_null<MassiveBody const*> secondary;
    Instant time;
    // The number of energy levels at which to draw equipotentials.
    int levels = 8;
    // The level of the equipotential going through the L₁ Lagrange point.
    int l1_level = 7;
    // Whether to draw an additional equipotential going through the L₂ point.
    bool show_l2_level = true;
  };

  struct Equipotentials {
    Equipotential<Barycentric, Navigation>::Lines lines;
    Parameters parameters;
  };

  GeometricPotentialPlotter(not_null<Ephemeris<Barycentric>*> ephemeris);

  // Cancel any computation in progress, causing the next call to
  // |RequestEquipotentials| to be processed as fast as possible.
  void Interrupt();

  // Sets the parameters that will be used for the computation of the next
  // equipotentials.
  void RequestEquipotentials(Parameters const& parameters);

  // The last value passed to |RequestEquipotentials|.
  std::optional<Parameters> const& last_parameters() const;

  // Sets |equipotentials()| to the latest computed equipotentials.
  void RefreshEquipotentials();

  Equipotentials const* equipotentials() const;

 private:

  absl::Status PlotEquipotentials(Parameters const& parameters);

  not_null<Ephemeris<Barycentric>*> const ephemeris_;

  std::optional<Parameters> last_parameters_;

  std::optional<Equipotentials> equipotentials_;

  mutable absl::Mutex lock_;
  jthread plotter_;

  // The |plotter_| is idle:
  // — if it is not joinable, e.g. because it was stopped by |Interrupt()|, or
  // — if it is done computing |next_equipotentials_| and has stopped or is
  //   about to stop executing.
  // If it is joined once idle (and joinable), it will not attempt to acquire
  // |lock_|.
  bool plotter_idle_ GUARDED_BY(lock_) = true;
  // |next_analysis_| is set by the |analyser_| thread; it is read and cleared
  // by the main thread.
  std::optional<Equipotentials> next_equipotentials_ GUARDED_BY(lock_);
  // |progress_of_next_analysis_| is set by the |analyser_| thread; it tracks
  // progress in computing |next_analysis_|.
  std::atomic<double> progress_of_next_analysis_ = 0;
};

}  // namespace internal

using internal::GeometricPotentialPlotter;

}  // namespace _geometric_potential_plotter
}  // namespace ksp_plugin
}  // namespace principia
