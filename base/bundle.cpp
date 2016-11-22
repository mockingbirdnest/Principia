#pragma once

#include <functional>
#include <queue>

#include "base/status.hpp"
#include "bundle.hpp"

namespace principia {
namespace base {

Bundle::Bundle(int threads) {
  {
    std::unique_lock<std::mutex> l(lock_);
    joining_ = false;
  }
  for (int i = 0; i < threads; ++i) {
    threads_.emplace_back([this]() {
      Task current;
      Status status;
      for (;;) {
        {
          std::unique_lock<std::mutex> l(lock_);
          if (!status.ok()) {
            status_ = status;
          }
          if (joining_ && tasks_.empty()) {
            return;
          }
          tasks_not_empty_or_joining_.wait(
              l,
              [this]() { return !tasks_.empty() || joining_; });
          if (tasks_.empty()) {
            return;
          }
          current = std::move(tasks_.front());
          tasks_.pop();
        }
        status = current();
      }
    });
  }
}

void Bundle::Add(Task task) {
  {
    std::unique_lock<std::mutex> l(lock_);
    tasks_.emplace(std::move(task));
  }
  tasks_not_empty_or_joining_.notify_one();
}

Status Bundle::Join() {
  {
    std::unique_lock<std::mutex> l(lock_);
    joining_ = true;
  }
  tasks_not_empty_or_joining_.notify_all();
  for (auto& thread : threads_) {
    thread.join();
  }
  std::unique_lock<std::mutex> l(lock_);
  return status_;
}

}  // namespace base
}  // namespace principia
