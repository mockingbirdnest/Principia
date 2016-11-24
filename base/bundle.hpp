#pragma once

#include <atomic>
#include <experimental/optional>
#include <functional>
#include <queue>
#include <memory>
#include <mutex>
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

  // If a |task| returns an erroneous |Status|, the |Bundle| is aborted and
  // |Join| returns that status.
  void Add(Task task);

 private:
  // The |workers_| |Toil|.  This function returns when
  // |(tasks_.empty() && terminate_on_empty_) || abort_|.
  void Toil();

  // The |workers_| have their |AbortRequested| set to |BundleShouldAbort|.
  bool BundleShouldAbort();

  // Sets |abort_| and notifies all on |tasks_not_empty_or_terminate_|.
  void Abort();

  std::mutex lock_;
  // Notified once when a task is available for execution.  Notified to all
  // waiting workers when either |terminate_on_empty_| or |abort_| is set.
  // Workers wait on this when no |tasks_| are available.
  std::condition_variable tasks_not_empty_or_terminate_;

  // If |abort_|, currently-running tasks should cooperatively abort, and
  // workers will terminate without considering queued tasks.
  std::atomic<bool> abort_;
  // Whether the workers should terminate when no tasks are available.  Set by
  // |Join|.
  bool terminate_on_empty_ GUARDED_BY(lock_);
  std::queue<Task> tasks_ GUARDED_BY(lock_);
  // The first erroneous status returned by a task, or |OK| if everything is
  // fine.
  Status status_ GUARDED_BY(lock_);
  std::vector<std::thread> workers_ GUARDED_BY(lock_);

  // A pointer to |AbortRequested| on the thread which created the |Bundle|.
  int const max_workers_;
  not_null<std::function<bool()>*> const master_abort_;
};

}  // namespace base
}  // namespace principia
