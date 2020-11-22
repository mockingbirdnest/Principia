#include "base/jthread.hpp"

#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

TEST(JThreadTest, StopToken) {
  auto f = [](stop_token st, int* const value, bool* const observed_stop) {
    while (!st.stop_requested()) {
      std::cout << ++(*value) << ' ' << std::flush;
      absl::SleepFor(absl::Milliseconds(10));
    }
    std::cout << std::endl;
    *observed_stop = true;
  };

  int i = 5;
  bool observed_stop = false;
  {
    jthread thread(f, &i, &observed_stop);
    absl::SleepFor(absl::Milliseconds(100));
  }
  EXPECT_TRUE(observed_stop);
  EXPECT_GT(i, 14);  // May flake.
}

}  // namespace base
}  // namespace principia