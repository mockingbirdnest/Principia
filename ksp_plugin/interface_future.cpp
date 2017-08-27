
#include "ksp_plugin/interface.hpp"

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
constexpr int window_size = 500;
int window_index = 0;
Time total_Δt;
Time min_Δt = Infinity<Time>();
Time max_Δt = -Infinity<Time>();
std::chrono::steady_clock::time_point time_of_last_promise_batch;
bool last_call_was_a_promise = false;
}  // namespace

std::future<void> const* principia__FutureCatchUpVessel(
    Plugin* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FutureCatchUpVessel> m({plugin, vessel_guid});
  if (!last_call_was_a_promise) {
    last_call_was_a_promise = true;
    time_of_last_promise_batch = std::chrono::steady_clock::now();
  }
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(vessel_guid);
  auto future = plugin->CatchUpVessel(vessel_guid);
  return m.Return(future.release());
}

void principia__FutureWait(std::future<void> const** const future) {
  journal::Method<journal::FutureWait> m({future}, {future});
  if (last_call_was_a_promise) {
    last_call_was_a_promise = false;
    auto const Δt = std::chrono::nanoseconds(std::chrono::steady_clock::now() -
                    time_of_last_promise_batch).count() * Nano(Second);
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
  TakeOwnership(future)->wait();
  return m.Return();
}

}  // namespace interface
}  // namespace principia
