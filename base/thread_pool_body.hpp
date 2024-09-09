#pragma once

#include "base/thread_pool.hpp"

#include <utility>

namespace principia {
namespace base {
namespace _thread_pool {
namespace internal {

// A helper function that is specialized for void because void is not really a
// type.
template<typename T>
void ExecuteAndSetValue(std::function<T()> const& function,
                        std::promise<T>& promise) {
  promise.set_value(function());
}

template<>
inline void ExecuteAndSetValue<void>(std::function<void()> const& function,
                                     std::promise<void>& promise) {
  function();
  promise.set_value();
}

template<typename T>
ThreadPool<T>::ThreadPool(std::int64_t const pool_size) {
  for (std::int64_t i = 0; i < pool_size; ++i) {
    threads_.emplace_back(std::bind(&ThreadPool::DequeueCallAndExecute, this));
  }
}

template<typename T>
ThreadPool<T>::~ThreadPool() {
  {
    absl::MutexLock l(&lock_);
    shutdown_ = true;
  }
  for (auto& thread : threads_) {
    thread.join();
  }
}

template<typename T>
std::future<T> ThreadPool<T>::Add(std::function<T()> function) {
  std::future<T> result;
  {
    absl::MutexLock l(&lock_);
    calls_.push_back({std::move(function), std::promise<T>()});
    result = calls_.back().promise.get_future();
  }
  return result;
}

template<typename T>
std::optional<std::future<T>>
ThreadPool<T>::TryAdd(std::function<T()> function) {
  std::optional<std::future<T>> result;

  // We use double locking to avoid contention when `TryAdd` fails.
  lock_.ReaderLock();
  std::int64_t const idle_threads = threads_.size() - busy_threads_;
  if (calls_.size() + 1 <= idle_threads) {
    lock_.ReaderUnlock();

    absl::MutexLock l(&lock_);
    // Queue the call iff we are sure that we have enough idle threads to be
    // able to schedule the call without blocking.
    std::int64_t const idle_threads = threads_.size() - busy_threads_;
    if (calls_.size() + 1 <= idle_threads) {
      calls_.push_back({std::move(function), std::promise<T>()});
      result = calls_.back().promise.get_future();
    }
  } else {
    lock_.ReaderUnlock();
  }

  return result;
}

template<typename T>
bool ThreadPool<T>::WaitUntilIdleFor(absl::Duration const duration) {
  absl::ReaderMutexLock l(&lock_);

  // Release this thread if there are enough idle threads to guarantee that the
  // call is able to schedule without blocking.
  auto const can_schedule_immediately = [this]() {
    std::int64_t const idle_threads = threads_.size() - busy_threads_;
    return calls_.size() + 1 <= idle_threads;
  };
  return lock_.AwaitWithTimeout(absl::Condition(&can_schedule_immediately),
                                duration);
}

template<typename T>
void ThreadPool<T>::DequeueCallAndExecute() {
  for (;;) {
    Call this_call;

    // Wait until either the queue contains an element or this class is shutting
    // down.
    {
      absl::MutexLock l(&lock_);

      auto const has_calls_or_shutdown = [this] {
        return shutdown_ || !calls_.empty();
      };
      lock_.Await(absl::Condition(&has_calls_or_shutdown));

      if (shutdown_) {
        break;
      }
      this_call = std::move(calls_.front());
      calls_.pop_front();
      ++busy_threads_;
    }

    // Execute the function without holding the `lock_` as it might take some
    // time.
    ExecuteAndSetValue(this_call.function, this_call.promise);

    {
      absl::MutexLock l(&lock_);
      --busy_threads_;
    }
  }
}

}  // namespace internal
}  // namespace _thread_pool
}  // namespace base
}  // namespace principia
