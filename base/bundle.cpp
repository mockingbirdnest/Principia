
#include "base/bundle.hpp"

#include <functional>
#include <list>

#include "base/map_util.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

void Bundle::Add(Task task) {
  absl::MutexLock l(&lock_);
  CHECK(!joining_);
  ++number_of_active_workers_;
  workers_.emplace_back(&Bundle::Toil, this, std::move(task));
}

Status Bundle::Join() {
  {
    absl::MutexLock l(&lock_);
    joining_ = true;
  }
  all_done_.WaitForNotification();
  JoinAll();
  absl::ReaderMutexLock status_lock(&status_lock_);
  return status_;
}

Status Bundle::JoinWithin(std::chrono::steady_clock::duration Δt) {
  {
    absl::MutexLock l(&lock_);
    joining_ = true;
  }
  if (!all_done_.WaitForNotificationWithTimeout(absl::FromChrono(Δt))) {
    absl::MutexLock l(&status_lock_);
    status_ = Status(Error::DEADLINE_EXCEEDED, "");
  }
  JoinAll();
  absl::ReaderMutexLock status_lock(&status_lock_);
  return status_;
}

Status Bundle::JoinBefore(std::chrono::system_clock::time_point t) {
  {
    absl::MutexLock l(&lock_);
    joining_ = true;
  }
  if (!all_done_.WaitForNotificationWithDeadline(absl::FromChrono(t))) {
    absl::MutexLock l(&status_lock_);
    status_ = Status(Error::DEADLINE_EXCEEDED, "");
  }
  JoinAll();
  absl::ReaderMutexLock status_lock(&status_lock_);
  return status_;
}

void Bundle::Toil(Task task) {
  Status const status = task();
  if (!status.ok()) {
    absl::MutexLock l(&status_lock_);
    status_.Update(status);
  }

  // Avoid locking when there are still active workers to reduce contention
  // during joining.
  if (--number_of_active_workers_ == 0) {
    absl::ReaderMutexLock l(&lock_);
    if (number_of_active_workers_ == 0 && joining_) {
      all_done_.Notify();
    }
  }
}

void Bundle::JoinAll() {
  absl::ReaderMutexLock l(&lock_);
  for (auto& worker : workers_) {
    worker.join();
  }
}

}  // namespace base
}  // namespace principia
