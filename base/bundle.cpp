#pragma once

#include "base/bundle.hpp"

#include <functional>

#include "base/map_util.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

thread_local std::function<bool()> AbortRequested = [] { return false; };

Bundle::TaskHandle::TaskHandle(std::shared_ptr<Task> task)
    : task_(std::move(task)) {}

Bundle::Bundle(int threads)
    : terminate_on_empty_(false),
      max_workers_(threads) {
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

Bundle::TaskHandle Bundle::Add(std::function<Status()> task_body) {
  if (workers_.size() < max_workers_) {
    workers_.emplace_back(&Bundle::Toil, this);
  }
  // Private constructor, cannot |make_shared|.
  TaskHandle result(std::make_shared<Task>(this, std::move(task_body)));
  {
    std::unique_lock<std::mutex> lock(lock_);
    CHECK(!terminate_on_empty_);
    tasks_.emplace_back(result.task_);
    auto end = tasks_.end();
    result.task_->set_queue_it(--end);
  }
  tasks_not_empty_or_terminate_.notify_one();
  return result;
}

void Bundle::Abort(TaskHandle const& handle) {
  not_null<Task*> task = handle.task_.get();
  task->Abort();
}

Status Bundle::JoinTask(TaskHandle const handle) {
  not_null<Task*> task = handle.task_.get();
  task->wait_for_termination();
  return task->status();
}

void Bundle::Toil() {
  std::shared_ptr<Task> current_task;
  for (;;) {
    {
      std::unique_lock<std::mutex> lock(lock_);
      if (terminate_on_empty_ && tasks_.empty()) {
        return;
      }
      tasks_not_empty_or_terminate_.wait(
          lock,
          [this] { return !tasks_.empty() || terminate_on_empty_; });
      if (tasks_.empty()) {
        return;
      }
      current_task = std::move(tasks_.front());
      tasks_.pop_front();
    }
    current_task->Execute();
    Status status = current_task->status();
    {
      std::unique_lock<std::mutex> lock(lock_);
      if (!status.ok() && status_.ok()) {
        status_ = status;
      }
    }
    current_task.reset();
  }
}

Bundle::Task::Task(not_null<Bundle*> bundle, std::function<Status()> body)
    : bundle_(bundle),
      body_(std::move(body)) {}

std::uint8_t Bundle::Task::execution_state() const {
  std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
  return execution_state_;
}

void Bundle::Task::set_queue_it(std::list<std::shared_ptr<Task>>::iterator it) {
  queue_it = it;
}

void Bundle::Task::Abort() {
  std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
  if (execution_state_ == Inactive) {
    bundle_->tasks_.erase(queue_it);
    status_ = Status(Error::ABORTED, "Task aborted in queue");
  }
  execution_state_ |= Aborted;
}

void Bundle::Task::Execute() {
  {
    std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
    CHECK_EQ(execution_state_, Inactive);
    execution_state_ |= Active;
    AbortRequested = [this] {
      // TODO(egg): timeout conditions.
      std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
      // the |!= 0| is there because of an MSVC warning (MSVC doesn't like
      // converting integers to booleans, even explicitly).  It also gives us
      // return type deduction.
      return (execution_state_ & Aborted) != 0;
    };
  }
  Status status = body_();
  {
    std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
    status_ = std::move(status);
    execution_state_ |= Terminated;
  }
  terminated_.notify_all();
}

Status Bundle::Task::status() const {
  std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
  CHECK_NE(execution_state_ & Terminated, 0) << execution_state_;
  return status_;
}

void Bundle::Task::wait_for_termination() {
  std::unique_lock<std::mutex> bundle_lock(bundle_->lock_);
  terminated_.wait(
      bundle_lock,
      [this] { return (execution_state_ & Terminated) != 0; });
}

}  // namespace base
}  // namespace principia
