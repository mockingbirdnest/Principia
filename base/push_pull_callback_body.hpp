#pragma once

#include "base/push_pull_callback.hpp"

#include <utility>

namespace principia {
namespace base {
namespace _push_pull_callback {
namespace internal {

template<typename Result, typename... Arguments>
bool PushPullCallback<Result, Arguments...>::Pull(Arguments&... arguments) {
  WaitUntilHasArgumentsOrShuttingDownAndLock();
  if (shutdown_) {
    return false;
  }
  std::tie(arguments...) = std::move(arguments_.value());
  arguments_.reset();
  lock_.Unlock();
  return true;
}

template<typename Result, typename... Arguments>
void PushPullCallback<Result, Arguments...>::Push(Result result) {
  absl::MutexLock l(&lock_);
  result_ = std::move(result);
}

template<typename Result, typename... Arguments>
std::function<Result(Arguments...)>
PushPullCallback<Result, Arguments...>::ToStdFunction() {
  return [this](Arguments const&... arguments) {
    Push(arguments...);
    return Pull();
  };
}

template<typename Result, typename ...Arguments>
void PushPullCallback<Result, Arguments...>::Shutdown() {
  absl::MutexLock l(&lock_);
  shutdown_ = true;
}

template<typename Result, typename... Arguments>
void PushPullCallback<Result, Arguments...>::Push(
    Arguments... arguments) {
  absl::MutexLock l(&lock_);
  arguments_ = std::tuple(std::move(arguments)...);
}

template<typename Result, typename... Arguments>
Result PushPullCallback<Result, Arguments...>::Pull() {
  WaitUntilHasResultAndLock();
  lock_.AssertHeld();
  Result result = result_.value();
  result_.reset();
  lock_.Unlock();
  return std::move(result);
}

template<typename Result, typename... Arguments>
void PushPullCallback<Result, Arguments...>::
    WaitUntilHasArgumentsOrShuttingDownAndLock() {
  auto has_arguments_or_shutting_down = [this]() {
    lock_.AssertReaderHeld();
    return shutdown_ || arguments_.has_value();
  };

  lock_.LockWhen(absl::Condition(&has_arguments_or_shutting_down));
}

template<typename Result, typename... Arguments>
void PushPullCallback<Result, Arguments...>::WaitUntilHasResultAndLock() {
  auto has_result = [this]() {
    lock_.AssertReaderHeld();
    return result_.has_value();
  };

  lock_.LockWhen(absl::Condition(&has_result));
}

template<typename T, typename Result, typename... Arguments>
PushPullExecutor<T, Result, Arguments...>::PushPullExecutor(Task task)
    : thread_([this, task = std::move(task)]() {
        auto const result = task(callback_.ToStdFunction());
        {
          absl::MutexLock l(&lock_);
          result_ = result;
        }
        callback_.Shutdown();
      }) {}

template<typename T, typename Result, typename... Arguments>
PushPullExecutor<T, Result, Arguments...>::~PushPullExecutor() {
  thread_.join();
}

template<typename T, typename Result, typename... Arguments>
PushPullCallback<Result, Arguments...>&
PushPullExecutor<T, Result, Arguments...>::callback() {
  absl::MutexLock l(&lock_);
  return callback_;
}

template<typename T, typename Result, typename... Arguments>
T PushPullExecutor<T, Result, Arguments...>::get() {
  absl::MutexLock l(&lock_);
  return std::move(result_.value());
}

}  // namespace internal
}  // namespace _push_pull_callback
}  // namespace base
}  // namespace principia
