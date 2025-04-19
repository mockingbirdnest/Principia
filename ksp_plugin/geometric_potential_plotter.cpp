#include "ksp_plugin/geometric_potential_plotter.hpp"

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

namespace principia {
namespace ksp_plugin {
namespace _geometric_potential_plotter {
namespace internal {

GeometricPotentialPlotter::GeometricPotentialPlotter(
    not_null<Ephemeris<Barycentric>*> const ephemeris)
    : ephemeris_(ephemeris) {}

void GeometricPotentialPlotter::Interrupt() {
  plotter_ = jthread();
  // We are single-threaded here, no need to lock.
  plotter_idle_ = true;
}

void GeometricPotentialPlotter::RequestEquipotentials(
    Parameters const& parameters) {
  last_parameters_ = parameters;
  absl::MutexLock l(&lock_);
  // Only process this request if there is no analysis in progress.
  if (plotter_idle_) {
    plotter_idle_ = false;
    plotter_ = MakeStoppableThread([this, parameters]() {
      auto const status = PlotEquipotentials(parameters);
      if (!status.ok() && !absl::IsCancelled(status)) {
        LOG(ERROR) << "Error while computing "
                   << parameters.primaries.front()->name()
                   << "-" << parameters.secondaries.front()->name()
                   << " equipotentials at " << parameters.time << ": "
                   << status;
      }
    });
  }
}

std::optional<GeometricPotentialPlotter::Parameters> const&
GeometricPotentialPlotter::last_parameters() const {
  return last_parameters_;
}

void GeometricPotentialPlotter::RefreshEquipotentials() {
  absl::MutexLock l(&lock_);
  if (next_equipotentials_.has_value()) {
    equipotentials_ = std::move(next_equipotentials_);
    next_equipotentials_.reset();
  }
}

GeometricPotentialPlotter::Equipotentials const*
GeometricPotentialPlotter::equipotentials() const {
  return equipotentials_.has_value() ? &*equipotentials_ : nullptr;
}

absl::Status GeometricPotentialPlotter::PlotEquipotentials(
    Parameters const& parameters) {
  auto result = LagrangeEquipotentials<Barycentric, Navigation>(ephemeris_)
                    .ComputeLines(parameters);

  absl::MutexLock l(&lock_);
  // We don’t reset `next_equipotentials_` unless `result.ok()`, so that if we
  // have a transient error, we keep the old ones until the problem goes away.
  if (result.ok()) {
    next_equipotentials_ = {{}, parameters};
    for (auto& [energy, lines] : result->lines) {
      next_equipotentials_->lines.insert(next_equipotentials_->lines.end(),
                                         std::make_move_iterator(lines.begin()),
                                         std::make_move_iterator(lines.end()));
    }
  }
  plotter_idle_ = true;
  return result.status();
}

}  // namespace internal
}  // namespace _geometric_potential_plotter
}  // namespace ksp_plugin
}  // namespace principia
