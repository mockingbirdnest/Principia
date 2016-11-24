#pragma once

#include <experimental/optional>
#include <functional>
#include <queue>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <vector>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

// Functions that can be cooperatively aborted should check |AbortRequested()|
// when appropriate.
extern thread_local std::function<bool()> AbortRequested;

// A thread pool with cooperative aborting (but no cooperative scheduling).
// We refer to the thread on which the |Bundle| is created as its master thread.
// The |Bundle| itself cooperatively aborts if |AbortRequested()| on its master
// thread.
class Bundle {
 public:
  using Task = std::function<Status()>;

  explicit Bundle(int workers);

  // Returns the first non-OK status encountered, or OK.  All worker threads are
  // joined; no calls to member functions may follow this call.
  Status Join();
  // Same as above, but aborts and returns |Error::DEADLINE_EXCEEDED| if it
  // fails to complete within the given interval.
  Status JoinWithin();
  // Same as above with absolute time.
  Status JoinBefore();

  // If a |task| returns an erroneous |Status|, the |Bundle| is aborted and
  // |Join| returns that status.
  void Add(Task task);

 private:
  // The |workers_| |Toil|.  This function returns when
  // |(tasks_.empty() && terminate_on_empty_) || abort_|.
  void Toil();

  // The |workers_| have their |AbortRequested| set to |BundleShouldAbort|.
  bool BundleShouldAbort();

  // If |status_| is erroneous, has no effect. Otherwise, sets |status_| to
  // |status| and notifies all on |tasks_not_empty_or_terminate_|.  |status|
  // should not be |OK|.
  void Abort(Status status);

  // Shared-locks |abort_lock_| and returns |!status_.ok()|.
  bool Aborting();

  // |abort_lock_| should not be held when locking |lock_|.
  std::mutex lock_;
  // Notified once when a task is available for execution.  Notified to all
  // waiting workers when either |terminate_on_empty_| or |abort_| is set.
  // Workers wait on this when no |tasks_| are available.
  std::condition_variable tasks_not_empty_or_terminate_;

  // Only lock in |Abort| and |Aborting|, and at the end of |Join|.
  std::shared_mutex abort_lock_;
  // If |!status_.ok()|, currently-running tasks should cooperatively abort, and
  // workers will terminate without considering queued tasks.  Set by |Abort|,
  // accessed by |Aborting|, returned by |Join|.
  Status status_ GUARDED_BY(abort_lock_);

  // Whether the workers should terminate when no tasks are available.  Set by
  // |Join|.
  bool terminate_on_empty_ GUARDED_BY(lock_);
  std::queue<Task> tasks_ GUARDED_BY(lock_);
  std::vector<std::thread> workers_ GUARDED_BY(lock_);

  // A pointer to |AbortRequested| on the thread which created the |Bundle|.
  int const max_workers_;
  not_null<std::function<bool()>*> const master_abort_;
};

}  // namespace base
}  // namespace principia
