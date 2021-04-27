#pragma once

#include "base/recurring_thread.hpp"

#include "gtest/gtest.h"

namespace principia {

using namespace std::chrono_literals;

namespace base {

class RecurringThreadTest : public ::testing::Test {
 protected:
  using ToyRecurringThread = RecurringThread<int, double>;

  static double PoolingGet(ToyRecurringThread& thread) {
    std::optional<double> output;
    do {
      output = thread.Get();
      std::this_thread::sleep_for(50us);
    } while (!output.has_value());
    return output.value();
  }
};

TEST_F(RecurringThreadTest, Basic) {
  auto add_one_half = [](int const input) {
    return static_cast<double>(input) + 0.5;
  };

  ToyRecurringThread thread(std::move(add_one_half), 1ms);
  thread.Start();

  thread.Put(3);
  {
    double const output = PoolingGet(thread);
    EXPECT_EQ(3.5, output);
  }

  EXPECT_FALSE(thread.Get().has_value());

  thread.Put(4);
  {
    double const output = PoolingGet(thread);
    EXPECT_EQ(4.5, output);
  }
}

}  // namespace base
}  // namespace principia
