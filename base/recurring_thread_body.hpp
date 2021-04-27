#pragma once

#include "base/recurring_thread.hpp"

#include <algorithm>

namespace principia {
namespace base {
namespace internal_recurring_thread {

void BaseRecurringThread::Start() {
  absl::MutexLock l(&jthread_lock_);
  if (!jthread_.joinable()) {
    jthread_ = MakeStoppableThread(
        [this]() { Status const status = RepeatedlyRunAction(); });
  }
}

void BaseRecurringThread::Stop() {
  absl::MutexLock l(&jthread_lock_);
  jthread_ = jthread();
}

BaseRecurringThread::BaseRecurringThread(std::chrono::milliseconds const period)
    : period_(period) {}

Status BaseRecurringThread::RepeatedlyRunAction() {
  for (std::chrono::steady_clock::time_point wakeup_time;;
       std::this_thread::sleep_until(wakeup_time)) {
    wakeup_time = std::chrono::steady_clock::now() + period_;
    RETURN_IF_STOPPED;

    RunAction();

    RETURN_IF_STOPPED;
  }
}

template<typename Input, typename Output>
RecurringThread<Input, Output>::RecurringThread(
    Action action,
    std::chrono::milliseconds const period)
    : BaseRecurringThread(period),
      action_(std::move(action)) {}

template<typename Input, typename Output>
void RecurringThread<Input, Output>::Put(Input input) {
  absl::MutexLock l(&input_output_lock_);
  input_ = std::move(input);
}

template<typename Input, typename Output>
std::optional<Output> RecurringThread<Input, Output>::Get() {
  absl::MutexLock l(&input_output_lock_);
  std::optional<Output> result;
  if (output_.has_value()) {
    std::swap(result, output_);
  }
  return result;
}

template<typename Input, typename Output>
Status RecurringThread<Input, Output>::RunAction() {
  std::optional<Input> input;
  {
    absl::MutexLock l(&input_output_lock_);
    if (!input_.has_value()) {
      // No input, let's wait for it to appear.
      return Status::OK;
    }
    std::swap(input, input_);
  }
  RETURN_IF_STOPPED;

  StatusOr<Output> status_or_output = action_(input.value());
  RETURN_IF_STOPPED;

  if (status_or_output.ok()) {
    absl::MutexLock l(&input_output_lock_);
    output_ = std::move(status_or_output.ValueOrDie());
  }

  return status_or_output.status();
}

template<typename Input>
RecurringThread<Input, void>::RecurringThread(
    Action action,
    std::chrono::milliseconds const period)
    : BaseRecurringThread(period),
      action_(std::move(action)) {}

template<typename Input>
void RecurringThread<Input, void>::Put(Input input) {
  absl::MutexLock l(&input_output_lock_);
  input_ = std::move(input);
}

template<typename Input>
Status RecurringThread<Input, void>::RunAction() {
  std::optional<Input> input;
  {
    absl::MutexLock l(&input_output_lock_);
    if (!input_.has_value()) {
      // No input, let's wait for it to appear.
      return Status::OK;
    }
    std::swap(input, input_);
  }
  RETURN_IF_STOPPED;

  return action_(input.value());
}

}  // namespace internal_recurring_thread
}  // namespace base
}  // namespace principia
