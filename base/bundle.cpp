#pragma once

#include "base/bundle.hpp"

#include <functional>
#include <list>

#include "base/map_util.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

thread_local std::function<bool()> AbortRequested = [] { return false; };

Bundle::Bundle(int threads)
    : abort_(false),
      terminate_on_empty_(false),
      max_workers_(threads),
      master_abort_(&AbortRequested) {
  workers_.reserve(max_workers_);
}

Status Bundle::Join() {
  {
    std::unique_lock<std::mutex> lock(lock_);
    CHECK(!terminate_on_empty_);
    terminate_on_empty_ = true;
  }
  tasks_not_empty_or_terminate_.notify_all();
  for (auto& thread : workers_) {
    thread.join();
  }
  std::unique_lock<std::mutex> lock(lock_);
  return status_;
}

void Bundle::Add(Task task) {
  {
    std::unique_lock<std::mutex> lock(lock_);
    CHECK(!terminate_on_empty_);
    if (workers_.size() < max_workers_) {
      workers_.emplace_back(&Bundle::Toil, this);
    }
    tasks_.emplace(task);
  }
  tasks_not_empty_or_terminate_.notify_one();
}

void Bundle::Toil() {
  Task current_task;
  AbortRequested = std::bind(&Bundle::BundleShouldAbort, this);
  for (;;) {
    {
      std::unique_lock<std::mutex> lock(lock_);
      tasks_not_empty_or_terminate_.wait(
          lock,
          [this] { return !tasks_.empty() || terminate_on_empty_ || abort_; });
      if (tasks_.empty() || abort_) {
        if (abort_) {
          LOG(ERROR) << std::this_thread::get_id() << " aborts";
        } else {
          LOG(ERROR) << std::this_thread::get_id() << " joins (no more tasks)";
        }
        return;
      }
      current_task = std::move(tasks_.front());
      tasks_.pop();
    }
    LOG(ERROR) << std::this_thread::get_id() << " executes a task";
    Status status = current_task();
    bool new_erroneous_status = false;
    {
      std::unique_lock<std::mutex> lock(lock_);
      if (!status.ok() && status_.ok()) {
        status_ = status;
        Abort();
      }
    }
  }
}

bool Bundle::BundleShouldAbort() {
  if (!abort_ && (*master_abort_)()) {
    Abort();
  }
  return abort_;
}

void Bundle::Abort() {
  abort_ = true;
  tasks_not_empty_or_terminate_.notify_all();
}

}  // namespace base
}  // namespace principia
