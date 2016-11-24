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
    : terminate_on_empty_(false),
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
  std::shared_lock<std::shared_mutex> abort_lock(abort_lock_);
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
          [this] {
            return !tasks_.empty() || terminate_on_empty_ || Aborting();
          });
      if (tasks_.empty() || Aborting()) {
        if (Aborting()) {
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
    if (!status.ok()) {
      Abort(status);
    }
  }
}

bool Bundle::BundleShouldAbort() {
  if (!Aborting() && (*master_abort_)()) {
    Abort(Status(Error::ABORTED, "abort requested on bundle master"));
  }
  return Aborting();
}

void Bundle::Abort(Status status) {
  std::unique_lock<std::shared_mutex> abort_lock(abort_lock_);
  CHECK(!status.ok());
  if (!status_.ok()) {
    // Already aborting.
    return;
  }
  status_ = status;
  tasks_not_empty_or_terminate_.notify_all();
}

bool Bundle::Aborting() {
  std::shared_lock<std::shared_mutex> abort_lock(abort_lock_);
  return !status_.ok();
}

}  // namespace base
}  // namespace principia
