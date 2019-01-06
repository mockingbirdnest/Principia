
#pragma once

#include <functional>
#include <list>
#include <optional>
#include <thread>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/synchronization/notification.h"
#include "absl/synchronization/mutex.h"
#include "base/status.hpp"

namespace principia {
namespace base {

// A bundle manages a number of threads that execute independently.  A thread is
// created for each call to |Add|.  When one of the |Join*| is called, no more
// calls to |Add| are allowed, and |Join*| returns the first error status (if
// any) produced by the tasks.
class Bundle final {
 public:
  using Task = std::function<Status()>;

  // If a |task| returns an erroneous |Status|, |Join| returns that status.
  void Add(Task task) LOCKS_EXCLUDED(lock_);

  // Returns the first non-OK status encountered, or OK.  All worker threads are
  // joined; no calls to member functions may follow this call.
  Status Join() LOCKS_EXCLUDED(lock_, status_lock_);
  // Same as above, but returns |Error::DEADLINE_EXCEEDED| if it fails to
  // complete within the given interval.
  Status JoinWithin(std::chrono::steady_clock::duration Δt)
      LOCKS_EXCLUDED(lock_, status_lock_);
  // Same as above with absolute time.
  Status JoinBefore(std::chrono::system_clock::time_point t)
      LOCKS_EXCLUDED(lock_, status_lock_);

 private:
  // Run on a separate thread to execute task and record its status.
  void Toil(Task const& task) LOCKS_EXCLUDED(lock_, status_lock_);

  void JoinAll() LOCKS_EXCLUDED(lock_);

  absl::Mutex status_lock_;
  Status status_ GUARDED_BY(status_lock_);

  absl::Mutex lock_;

  // Whether |Join| has been called.  When set to true, |Add| should not be
  // called.
  std::atomic_bool joining_ = false;
  absl::Notification all_done_;

  // The number of workers currently executing.  Can only be incremented when
  // |joining_| is false.
  std::atomic_int number_of_active_workers_ = 0;
  std::list<std::thread> workers_ GUARDED_BY(lock_);

  static_assert(std::atomic_bool::is_always_lock_free, "bool not lock-free");
  static_assert(std::atomic_int::is_always_lock_free, "int not lock-free");
};

}  // namespace base
}  // namespace principia
