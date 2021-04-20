#pragma once

#include "base/recurring_thread.hpp"

namespace principia {
namespace base {
namespace internal_recurring_thread {

template<typename Input, typename Output>
RecurringThread<Input, Output>::RecurringThread(
    Action action,
    std::chrono::milliseconds const period)
    : action_(std::move(action)),
      period_(period),
      jthread_(MakeStoppableThread(
          [this]() { Status const status = RepeatedlyRunAction(); })) {}

template<typename Input, typename Output>
void RecurringThread<Input, Output>::Put(Input input) {
  absl::MutexLock l(&lock_);
  input_ = std::move(input);
}

template<typename Input, typename Output>
std::optional<Output> RecurringThread<Input, Output>::Get() {
  absl::MutexLock l(&lock_);
  std::optional<Output> result;
  if (output_.has_value()) {
    std::swap(result, output_);
  }
  return result;
}

template<typename Input, typename Output>
Status RecurringThread<Input, Output>::RepeatedlyRunAction() {
  for (std::chrono::steady_clock::time_point wakeup_time;;
       std::this_thread::sleep_until(wakeup_time)) {
    wakeup_time = std::chrono::steady_clock::now() + period_;
    RETURN_IF_STOPPED;

    std::optional<Input> input;
    {
      absl::MutexLock l(&lock_);
      if (!input_.has_value()) {
        // No input, let's wait for it to appear.
        continue;
      }
      std::swap(input, input_);
    }
    RETURN_IF_STOPPED;

    Output output = action_(input.value());
    RETURN_IF_STOPPED;

    {
      absl::MutexLock l(&lock_);
      output_ = std::move(output);
    }
    RETURN_IF_STOPPED;
  }
}

}  // namespace internal_recurring_thread
}  // namespace base
}  // namespace principia
