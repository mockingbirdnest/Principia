
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
  joining_ = true;
  all_done_.WaitForNotification();
  JoinAll();
  absl::ReaderMutexLock status_lock(&status_lock_);
  return status_;
}

Status Bundle::JoinWithin(std::chrono::steady_clock::duration Δt) {
  joining_ = true;
  if (!all_done_.WaitForNotificationWithTimeout(absl::FromChrono(Δt))) {
    absl::MutexLock l(&status_lock_);
    status_ = Status(Error::DEADLINE_EXCEEDED, "bundle deadline exceeded");
  }
  JoinAll();
  absl::ReaderMutexLock status_lock(&status_lock_);
  return status_;
}

Status Bundle::JoinBefore(std::chrono::system_clock::time_point t) {
  joining_ = true;
  if (!all_done_.WaitForNotificationWithDeadline(absl::FromChrono(t))) {
    absl::MutexLock l(&status_lock_);
    status_ = Status(Error::DEADLINE_EXCEEDED, "bundle deadline exceeded");
  }
  JoinAll();
  absl::ReaderMutexLock status_lock(&status_lock_);
  return status_;
}

void Bundle::Toil(Task const& task) {
  Status const status = task();

  // Avoid locking if the task succeeded: it cannot affect the overall status.
  if (!status.ok()) {
    absl::MutexLock l(&status_lock_);
    status_.Update(status);
  }

  // No locking unless this is the last task and we are joining.  This
  // avoids contention during joining.  Note that if |joining_| is true we know
  // that |number_of_active_workers_| cannot increase.
  --number_of_active_workers_;
  if (joining_ && number_of_active_workers_ == 0) {
    all_done_.Notify();
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
