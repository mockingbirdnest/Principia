#pragma once

#include "base/reanimator.hpp"

#include <memory>
#include <utility>

#include "absl/log/check.h"
#include "absl/log/die_if_null.h"
#include "base/stoppable_thread.hpp"

namespace principia {
namespace base {
namespace _reanimator {
namespace internal {

using namespace principia::base::_stoppable_thread;

template<typename Key>
Reanimator<Key>::Reanimator(Action action) : action_(std::move(action)) {}

template<typename Key>
void Reanimator<Key>::Start() {
  absl::MutexLock l(&jthread_lock_);
  if (!jthread_.joinable()) {
    jthread_ = MakeStoppableThread([this]() { RepeatedRunActions(); });
  }
}

template<typename Key>
void Reanimator<Key>::Stop() {
  absl::MutexLock l(&jthread_lock_);

  // Force the thread to exit once all the queued actions have been executed.
  // Also ensure that after this point it's not possible to enqueue new actions.
  {
    absl::MutexLock l(&lock_);
    stopped_ = true;
    jthread_must_exit_ = true;
  }

  // Mark the thread as stopped and wait for all the queued actions to complete.
  // This includes the best-effort ones.  The stop request is so that actions
  // that have a `RETURN_IF_STOPPED` can finish quickly.
  jthread_.request_stop();
  jthread_.join();
}

template<typename Key>
void Reanimator<Key>::RunBestEffort(Key const& key) {
  absl::MutexLock l(&lock_);
  CHECK(!stopped_);

  // Queue the call.  The `PendingRun` is just owned by the queue.
  queue_.emplace(key, std::make_shared<PendingRun>());
}

template<typename Key>
typename Reanimator<Key>::Handle Reanimator<Key>::RunGuaranteed(
    Key const& key) {
  absl::MutexLock l(&lock_);
  CHECK(!stopped_);

  // Queue the call.  The `Notification` will be notified once the action has
  // run.
  auto const handle = std::make_shared<PendingRun>(
      PendingRun{.done = std::make_unique<absl::Notification>()});
  queue_.emplace(key, handle);

  // The `PendingRun` is co-owned by the queue and the caller of this function.
  return handle;
}

template<typename Key>
void Reanimator<Key>::Cancel(Key const& before_key) {
  absl::MutexLock l(&jthread_lock_);
  {
    absl::MutexLock l(&lock_);
    for (auto it = queue_.begin(); it != queue_.end();) {
      auto const& [key, run] = *it;
      if (run->done != nullptr) {
        // A guaranteed run.  Don't touch it, and don't touch any of the
        // following runs.
        return;
      } else if (key < before_key) {
        it = queue_.erase(it);
        // If the queue is now empty, we will mark the thread as stopped so
        // it should exit as soon as the current action is done.
        if (it == queue_.end()) {
          jthread_must_exit_ = true;
          break;
        }
      } else {
        // Here `key >= before_key`, so we are done.
        return;
      }
    }
  }

  // Mark the thread as stopped so that, if there is a running action with a
  // `RETURN_IF_STOPPED`, it can complete quickly.
  jthread_.request_stop();
  jthread_.join();
  {
    absl::MutexLock l(&lock_);
    jthread_must_exit_ = false;
  }

  // The thread that we stopped is now gone, create a new one.
  jthread_ = MakeStoppableThread([this]() { RepeatedRunActions(); });
}

template<typename Key>
absl::Status Reanimator<Key>::Wait(Handle const handle,
                                   ProgressCallback progress_callback) {
  // This object won't go away since we hold `handle`.
  auto const& pending_run = *ABSL_DIE_IF_NULL(handle);

  // Insert my progress callback, if any, and retain an iterator on it.
  std::optional<ProgressCallbacks::const_iterator> my_progress_callback;
  if (progress_callback != nullptr) {
    absl::MutexLock l(&lock_);
    progress_callbacks_.push_front(std::move(progress_callback));
    my_progress_callback = progress_callbacks_.begin();
  }

  ABSL_DIE_IF_NULL(pending_run.done)->WaitForNotification();

  // Remove my progress callback.
  if (my_progress_callback.has_value()) {
    absl::MutexLock l(&lock_);
    progress_callbacks_.erase(my_progress_callback.value());
  }

  return pending_run.status;
}

template<typename Key>
void Reanimator<Key>::RepeatedRunActions() {
  auto queue_not_empty_or_must_exit = [this] {
    lock_.AssertReaderHeld();
    return !queue_.empty() || jthread_must_exit_;
  };

  absl::MutexLock l(&lock_);
  for (;;) {
    lock_.Await(absl::Condition(&queue_not_empty_or_must_exit));

    if (queue_.empty()) {
      CHECK(jthread_must_exit_);
      break;
    }

    // Run the action at the end of the queue (the one with the greatest key).
    // Note that the run is still in the queue when the action is executed.  It
    // is important to retain a copy of `handle` because the queue might change
    // (due to `Cancel`) while we run the action.
    auto const [key, handle] = *queue_.rbegin();

    lock_.Unlock();
    handle->status = action_(key);
    lock_.Lock();

    // Run the progress callbacks.  This happens before unblocking the waiters
    // so that the waiters see their own progress.
    for (auto const& progress_callback : progress_callbacks_) {
      progress_callback(key, handle->status);
    }

    // Unblock the waiters.
    if (handle->done != nullptr) {
      handle->done->Notify();
    }

    // Remove the pending run from the queue.  We cannot hold an iterator across
    // the run of the action because the `queue_` could be modified while we
    // don't hold the `lock_`, including deleting the pending run that we are
    // executing.  Instead, we look explicitly for our key and handle.
    auto const [it1, it2] = queue_.equal_range(key);
    for (auto it = it1; it != it2; ++it) {
      if (it->second == handle) {
        // The `PendingRun` will be destroyed once all its handles are gone.
        queue_.erase(it);
        break;
      }
    }
  }
}

}  // namespace internal
}  // namespace _reanimator
}  // namespace base
}  // namespace principia
