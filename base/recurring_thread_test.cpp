#include "base/recurring_thread.hpp"

#include "gtest/gtest.h"

namespace principia {

using namespace std::chrono_literals;

namespace base {

using namespace principia::base::_recurring_thread;

class RecurringThreadTest : public ::testing::Test {
 protected:
  using ToyRecurringThread1 = RecurringThread<int>;
  using ToyRecurringThread2 = RecurringThread<int, double>;

  static double PollingGet(ToyRecurringThread2& thread) {
    std::optional<double> output;
    do {
      output = thread.Get();
      std::this_thread::sleep_for(50us);
    } while (!output.has_value());
    return output.value();
  }
};

TEST_F(RecurringThreadTest, Result) {
  auto add_one_half = [](int const input) {
    return static_cast<double>(input) + 0.5;
  };

  ToyRecurringThread2 thread(std::move(add_one_half), 1ms);
  thread.Start();

  thread.Put(3);
  {
    double const output = PollingGet(thread);
    EXPECT_EQ(3.5, output);
  }

  EXPECT_FALSE(thread.Get().has_value());

  thread.Put(4);
  {
    double const output = PollingGet(thread);
    EXPECT_EQ(4.5, output);
  }
}

TEST_F(RecurringThreadTest, NoResult) {
  double value = 0.0;

  auto add_one_half = [&value](int const input) {
    value = static_cast<double>(input) + 0.5;
    return absl::OkStatus();
  };

  ToyRecurringThread1 thread(std::move(add_one_half), 1ms);
  thread.Start();

  thread.Put(3);
  do {
    std::this_thread::sleep_for(50us);
  } while (value != 3.5);
}

}  // namespace base
}  // namespace principia
