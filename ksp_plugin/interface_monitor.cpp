
#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <string>
#include <type_traits>

#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using quantities::Infinity;
using quantities::Time;
using quantities::si::Nano;
using quantities::si::Second;

namespace {
constexpr int monitor_count = 64;
constexpr int window_size = 500;

struct Monitor {
  std::string* name = nullptr;
  std::chrono::steady_clock::time_point start_time{};
  bool is_running = false;

  int window_index = 0;
  Time total_Δt;
  Time min_Δt = Infinity<Time>;
  Time max_Δt = -Infinity<Time>;
};

static_assert(
    std::is_trivially_destructible<std::array<Monitor, monitor_count>>::value,
    "An array of |Monitor|s should be trivially destructible");
std::array<Monitor, monitor_count> monitors{};
}  // namespace

// No journalling to avoid measuring overhead from that; these functions have no
// side effects aside from logging and changing the internal state of the
// monitors.

void __cdecl principia__MonitorSetName(int const i, char const* const name) {
  Monitor& monitor = monitors[i];
  if (monitor.name == nullptr) {
    monitor.name = new std::string(name);
  }
}

void __cdecl principia__MonitorStart(int const i) {
  Monitor& monitor = monitors[i];
  if (!monitor.is_running) {
    monitor.is_running = true;
    monitor.start_time = std::chrono::steady_clock::now();
  }
}

void __cdecl principia__MonitorStop(int const i) {
  Monitor& monitor = monitors[i];
  if (monitor.is_running) {
    monitor.is_running = false;
    auto const Δt =
        std::chrono::nanoseconds(
            std::chrono::steady_clock::now() - monitor.start_time).count() *
        Nano(Second);
    monitor.min_Δt = std::min(monitor.min_Δt, Δt);
    monitor.max_Δt = std::max(monitor.max_Δt, Δt);
    monitor.total_Δt += Δt;
    ++monitor.window_index %= window_size;
    if (monitor.window_index == 0) {
      LOG(INFO) << "[Monitor " << i
                << (monitor.name == nullptr ? "" : (": " + *monitor.name))
                << "] min = " << monitor.min_Δt << ", max = " << monitor.max_Δt
                << u8", μ = " << monitor.total_Δt / window_size;
      monitor.min_Δt = Infinity<Time>;
      monitor.max_Δt = -Infinity<Time>;
      monitor.total_Δt = Time();
    }
  }
}

}  // namespace interface
}  // namespace principia
