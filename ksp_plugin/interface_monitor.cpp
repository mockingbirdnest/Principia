
#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <chrono>

#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using quantities::Infinity;
using quantities::Time;
using quantities::si::Nano;
using quantities::si::Second;

namespace {
std::chrono::steady_clock::time_point monitor_start_time;
bool is_monitor_running = false;
}  // namespace

// No journalling to avoid overhead from that; these functions have no side
// effects aside from logging changing and the internal state of the monitor.

void principia__MonitorStart() {
  if (!is_monitor_running) {
    is_monitor_running = true;
    monitor_start_time = std::chrono::steady_clock::now();
  }
}

void principia__MonitorStop() {
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

}  // namespace interface
}  // namespace principia
