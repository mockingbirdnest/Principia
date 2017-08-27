
#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <chrono>

#include "journal/method.hpp"
#include "journal/profiles.hpp"

namespace principia {
namespace interface {

using quantities::Infinity;
using quantities::Time;
using quantities::si::Nano;
using quantities::si::Second;

namespace {
std::chrono::steady_clock::time_point monitor_start_time;
bool is_monitor_running = false;

void StartMonitoring() {
  if (!is_monitor_running) {
    is_monitor_running = true;
    monitor_start_time = std::chrono::steady_clock::now();
  }
}

void StopMonitoring() {
  constexpr int window_size = 500;
  static int window_index = 0;
  static Time total_Δt;
  static Time min_Δt = Infinity<Time>();
  static Time max_Δt = -Infinity<Time>();
  if (is_monitor_running) {
    is_monitor_running = false;
    auto const Δt =
        std::chrono::nanoseconds(
            std::chrono::steady_clock::now() - monitor_start_time).count() *
        Nano(Second);
    min_Δt = std::min(min_Δt, Δt);
    max_Δt = std::max(max_Δt, Δt);
    total_Δt += Δt;
    ++window_index %= window_size;
    if (window_index == 0) {
      LOG(INFO) << "min = " << min_Δt << ", max = " << max_Δt
                << u8", μ = " << total_Δt / window_size;
      min_Δt = Infinity<Time>();
      max_Δt = -Infinity<Time>();
      total_Δt = Time();
    }
  }
}

}  // namespace

std::future<void> const* principia__FutureCatchUpVessel(
    Plugin* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FutureCatchUpVessel> m({plugin, vessel_guid});
  void StartMonitoring();
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(vessel_guid);
  auto future = plugin->CatchUpVessel(vessel_guid);
  return m.Return(future.release());
}

void principia__FutureWait(std::future<void> const** const future) {
  journal::Method<journal::FutureWait> m({future}, {future});
  StopMonitoring();
  TakeOwnership(future)->wait();
  return m.Return();
}

}  // namespace interface
}  // namespace principia
