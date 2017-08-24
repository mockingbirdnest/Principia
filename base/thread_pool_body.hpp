
#pragma once

#include "base/thread_pool.hpp"

namespace principia {
namespace base {

template<typename T>
ThreadPool<T>::ThreadPool(std::int64_t const pool_size) {
  for (std::int64_t i = 0; i < pool_size; ++i) {
    threads_.emplace_back(std::bind(&DequeCallAndExecute, this));
  }
}

template<typename T>
ThreadPool<T>::~ThreadPool() {
  {
    std::lock_guard<std::mutex> l(&lock_);
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
    std::lock_guard<std::mutex> l(&lock_);
    calls_.emplace_back({std::move(function), std::promise<T>()});
    result = calls_.back().promise.get_future();
  }
  has_calls_.notify_one();
  return result;
}

template<typename T>
void ThreadPool<T>::DequeueCallAndExecute() {
  for (;;) {
    Call this_call;
    {
      std::unique_lock<std::mutex> l(&lock_);
      has_calls_.wait(l, [] { return shutdown_ || !calls_.empty(); });
      if (shutdown_) {
        break;
      }
      this_call = calls_.pop_front();
    }
    this_call.promise.set_value(this_call.function());
  }
}

}  // namespace base
}  // namespace principia
