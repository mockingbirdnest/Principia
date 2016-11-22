#pragma once

#include "base/bundle.hpp"

#include <functional>
#include <queue>

#include "base/status.hpp"

namespace principia {
namespace base {

Bundle::Bundle(int threads) {
  {
    std::unique_lock<std::mutex> l(lock_);
    joining_ = false;
  }
  for (int i = 0; i < threads; ++i) {
    threads_.emplace_back([this]() {
      std::shared_ptr<Task> current_task = nullptr;
      for (;;) {
        {
          std::unique_lock<std::mutex> l(lock_);
          if (joining_ && tasks_.empty()) {
            return;
          }
          tasks_not_empty_or_joining_.wait(
              l,
              [this]() { return !tasks_.empty() || joining_; });
          if (tasks_.empty()) {
            return;
          }
          current_task = std::move(tasks_.front());
          tasks_.pop_front();
          current_task->in_queue_ = std::experimental::nullopt;
          current_task->abort_requested_ = &AbortRequested;
        }
        Status const status = current_task->body_();
        {
          std::unique_lock<std::mutex> l(lock_);
          current_task->abort_requested_ = nullptr;
          current_task->status_ = status;
          if (!status.ok() && status_.ok()) {
            status_ = status;
          }
        }
        current_task->done_.notify_all();
        current_task = nullptr;
      }
    });
  }
}

std::shared_ptr<Bundle::Task> Bundle::Add(std::function<Status()> task) {
  auto result = std::make_shared<Task>(task);
  {
    std::unique_lock<std::mutex> l(lock_);
    CHECK(!joining_);
    tasks_.emplace_back(result);
    tasks_.back();
  }
  tasks_not_empty_or_joining_.notify_one();
  return result;
}

Status Bundle::Join() {
  {
    std::unique_lock<std::mutex> l(lock_);
    CHECK(!joining_);
    joining_ = true;
  }
  tasks_not_empty_or_joining_.notify_all();
  for (auto& thread : threads_) {
    thread.join();
  }
  std::unique_lock<std::mutex> l(lock_);
  return status_;
}

void Bundle::Abort(not_null<Task*> task) {
  std::unique_lock<std::mutex> l(lock_);
  if (task->in_queue_) {
    tasks_.erase(*task->in_queue_);
    task->in_queue_ = std::experimental::nullopt;
  } else if (task->abort_requested_ != nullptr) {
    *task->abort_requested_ = true;
  }
  // TODO(egg): Should we do something, e.g., set task->status_ to aborted if it
  // terminated normally then was aborted?
}

Status Bundle::JoinTask(not_null<Task*> task) {
  std::unique_lock<std::mutex> l(lock_);
  // TODO(egg): this should use has_value rather than static_cast once we have
  // a better optional.
  task->done_.wait(l, [task]() { return static_cast<bool>(task->status_); });
  return *task->status_;
}

Bundle::Task::Task(std::function<Status()> body) : body_(std::move(body)) {}

}  // namespace base
}  // namespace principia
