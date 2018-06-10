
#pragma once

#include <functional>
#include <queue>
#include <memory>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include <thread>
#include <vector>

#include "base/macros.hpp"
#include "base/monostable.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

// Functions that can be cooperatively aborted should check |AbortRequested()|
// when appropriate.  The function is thread-safe.
extern thread_local std::function<bool()> AbortRequested;

// A thread pool with cooperative aborting (but no cooperative scheduling).
// We refer to the thread on which the |Bundle| is created as its master thread.
// The |Bundle| itself cooperatively aborts if |AbortRequested()| on its master
// thread.
class Bundle final {
 public:
  using Task = std::function<Status()>;

  explicit Bundle(int workers);

  // Returns the first non-OK status encountered, or OK.  All worker threads are
  // joined; no calls to member functions may follow this call.
  Status Join();
  // Same as above, but aborts and returns |Error::DEADLINE_EXCEEDED| if it
  // fails to complete within the given interval.
  Status JoinWithin(std::chrono::steady_clock::duration Δt);
  // Same as above with absolute time.
  Status JoinBefore(std::chrono::steady_clock::time_point t);

  // If a |task| returns an erroneous |Status|, the |Bundle| is aborted and
  // |Join| returns that status.
  void Add(Task task);

 private:
  // The |workers_| |Toil|.  This function returns when
  // |(tasks_.empty() && !wait_on_empty_) ||Aborting()|.
  void Toil();

  // The |workers_| have their |AbortRequested| set to |BundleShouldAbort|.
  bool BundleShouldAbort();

  // If |status_| is erroneous, has no effect. Otherwise, sets |status_| to
  // |status| and notifies all on |tasks_not_empty_or_terminate_|.  |status|
  // should not be |OK|.
  void Abort(Status status);

  // Thread-safe |!status_.ok()|.
  bool Aborting();

  // Thread-safe |deadline_ && std::chrono::steady_clock::now() > deadline_|.
  bool DeadlineExceeded();

  // |status_lock_| should not be held when locking |lock_|.
  std::mutex lock_;
  // Notified once when a task is available for execution.  Notified to all
  // waiting workers when either |wait_on_empty_| is flopped or |status_| is set
  // to an error.
  // Workers wait on this when no |tasks_| are available.
  std::condition_variable tasks_not_empty_or_terminate_;
  std::shared_mutex status_lock_;
  // If |!status_.ok()|, currently-running tasks should cooperatively abort, and
  // workers will terminate without considering queued tasks.  Set by |Abort|,
  // accessed by |Aborting|, returned by |Join|.
  Status status_ GUARDED_BY(status_lock_);
  std::optional<std::chrono::steady_clock::time_point> deadline_
      GUARDED_BY(status_lock_);

  // Whether the workers should terminate when no tasks are available.  Set by
  // |Join|.
  Monostable wait_on_empty_ GUARDED_BY(lock_);

  std::queue<Task> tasks_ GUARDED_BY(lock_);
  std::vector<std::thread> workers_ GUARDED_BY(lock_);

  int const max_workers_;
  // A pointer to |AbortRequested| on the thread which created the |Bundle|.
  not_null<std::function<bool()>*> const master_abort_;
};

}  // namespace base
}  // namespace principia
