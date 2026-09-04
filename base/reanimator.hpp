#pragma once

#include <functional>
#include <list>
#include <memory>
#include <thread>

#include "absl/base/thread_annotations.h"
#include "absl/container/btree_map.h"
#include "absl/synchronization/notification.h"
#include "absl/synchronization/mutex.h"
#include "absl/status/status.h"

namespace principia {
namespace base {
namespace _reanimator {
namespace internal {

// A helper class for reanimation.  It is essentially a `jthread` with a queue
// ordered by `Key`.  Actions may be executed best-effort (in which case they
// may be cancelled by `Cancel`) or with guaranteed execution (in which case
// they cannot be cancelled and a handle is returned that makes it possible to
// wait for the action to complete).
template<typename Key>
class Reanimator {
  struct PendingRun;
  using Handle = std::shared_ptr<PendingRun>;

 public:
  using Action = std::function<absl::Status(Key const&)>;
  using ProgressCallback =
      std::function<void(Key const& key, absl::Status const& status)>;

  explicit Reanimator(Action action);

  // Creates the internal thread which may start running queued actions.
  // Idempotent.
  void Start();

  // Waits for all queued actions (including the best-effort ones) to complete
  // and destroys the internal thread.  No action may be queued after this call.
  void Stop();

  // Queues a run of the action with the given `key`.  This run may never happen
  // if it is cancelled.
  void RunBestEffort(Key const& key);

  // Queues a run of the action with the given `key`.  This run is sure to
  // happen and the caller may use the returned handle to wait for its
  // completion.
  Handle RunGuaranteed(Key const& key);

  // Cancels all the best-effort runs with a key strictly less than
  // `before_key`.  This may include the action being executed.
  void Cancel(Key const& before_key);

  // Waits for the run with the given `handle` to complete and returns its
  // status.  The `progress_callback`, if any, is executed each time an action
  // completes.
  absl::Status Wait(Handle handle,
                    ProgressCallback progress_callback = nullptr);

 private:
  // This is a list because we retain iterators in it, so we need iterator
  // stability.
  using ProgressCallbacks = std::list<ProgressCallback>;

  // The description of a queued or executing run.  `done` is null iff the
  // run is best-effort.
  struct PendingRun {
    std::unique_ptr<absl::Notification> done;
    absl::Status status;
  };

  // The execution loop, run on the `jthread_`.
  void RepeatedRunActions();

  // Construction parameter.
  Action const action_;

  absl::Mutex queue_lock_;

  // Set to `true` when `Stop` has been called.
  bool stopping_ = false ABSL_GUARDED_BY(queue_lock_);

  // The progress callbacks of all the waiters.  They are executed each time an
  // action completes.
  ProgressCallbacks progress_callbacks_ ABSL_GUARDED_BY(queue_lock_);

  // The queue is ordered by increasing keys, and the first run to execute is
  // the one with the greatest key.  The run being executed is still present in
  // the queue.
  absl::btree_multimap<Key, Handle> queue_ ABSL_GUARDED_BY(queue_lock_);

  absl::Mutex jthread_lock_;
  std::jthread jthread_ ABSL_GUARDED_BY(queue_lock_);
};

}  // namespace internal

using internal::Reanimator;

}  // namespace _reanimator
}  // namespace base
}  // namespace principia

#include "base/reanimator_body.hpp"