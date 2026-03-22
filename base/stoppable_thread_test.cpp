#include "base/stoppable_thread.hpp"

#include <chrono>
#include <ostream>
#include <thread>

#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using namespace principia::base::_stoppable_thread;
using namespace std::chrono_literals;

TEST(JThreadTest, ThisJThread) {
  bool observed_stop = false;

  auto sleepy_worker = MakeStoppableThread(
      [](bool* const observed_stop) {
        for (;;) {
          absl::SleepFor(absl::Milliseconds(10));
          if (this_stoppable_thread::get_stop_token().stop_requested()) {
            std::cout << "Sleepy worker is requested to stop\n";
            *observed_stop = true;
            return;
          }
          std::cout << "Sleepy worker goes back to sleep\n";
        }
      },
      &observed_stop);

  std::this_thread::sleep_for(30ms);
  std::cout << "Requesting stop of sleepy worker\n";
  sleepy_worker.request_stop();
  sleepy_worker.join();
  EXPECT_TRUE(observed_stop);
}

}  // namespace base
}  // namespace principia
