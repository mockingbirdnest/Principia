#pragma once

#include "absl/log/die_if_null.h"
#include "base/reanimator.hpp"
#include "base/stoppable_thread.hpp"

namespace principia {
namespace base {
namespace _reanimator {
namespace internal {

template<typename Key>
Reanimator<Key>::Reanimator(Action action) : action_(std::move(action)) {}

template<typename Key>
void Reanimator<Key>::Start() {
  absl::MutexLock l(&jthread_lock_);
  if (!jthread_.joinable()) {
    jthread_ = MakeStoppableThread([this]() { Loop(); });
  }
}

template<typename Key>
void Reanimator<Key>::Stop() {
  auto queue_empty = [this] {
    queue_lock_.AssertReaderHeld();
    return queue_.empty();
  };

  // Wait for all the calls to complete.  This includes the asynchronous calls.
  // Note that after this point it's not possible to enqueue new runs.
  absl::MutexLock l(&queue_lock_);
  stopping_ = true;
  queue_lock_.Await(absl::Condition(&queue_empty));

  absl::MutexLock l(&jthread_lock_);
  jthread_ = std::jthread();
}

template<typename Key>
void Reanimator<Key>::RunAsynchronously(Key const& key) {
  absl::MutexLock l(&queue_lock_);
  CHECK(!stopping_);
  queue_.emplace(key, PendingRun{});
}

template<typename Key>
typename Reanimator<Key>::Handle Reanimator<Key>::RunSynchronously(
    Key const& key) {
  absl::MutexLock l(&queue_lock_);
  CHECK(!stopping_);

  // Queue the call.  The `Notification` will be notified once it has run.
  auto const handle =
      std::make_shared<PendingRun>({.done = std::make_unique<Notification>()});
  queue_.emplace(key, handle);

  // The `PendingRun` is co-owned by the queue and the caller of this function.
  return handle;
}

template<typename Key>
void Reanimator<Key>::Cancel(Key const& before_key) {
  absl::MutexLock l(&queue_lock_);
  for (auto it = queue_.begin(); it != queue_.end();) {
    auto const& [key, run] = *it;
    if (run.done != nullptr) {
      return;
    } else if (key < before_key) {
      it = queue_.erase(it);
      if (it == queue_.end()) {
        // If the queue is now empty, kill the running thread.  Remember that
        // the action being executed is still in the queue.
        absl::MutexLock l(&jthread_lock_);
        jthread_ = std::jthread();
      }
    } else {
      return;
    }
  }
}

template<typename Key>
absl::Status Reanimator<Key>::Wait(Handle const handle,
                                   ProgressCallback progress_callback) {
  // This object won't go away since we hold `handle`.
  auto const& pending_run = *ABSL_DIE_IF_NULL(handle);

  // Insert our progress callback, if any, and retain an iterator on it.
  std::optional<ProgressCallback::const_iterator> my_progress_callback;
  if (progress_callback != nullptr) {
    absl::MutexLock l(&queue_lock_);
    progress_callbacks_.push_front(std::move(progress_callback));
    my_progress_callback = progress_callbacks_.begin();
  }

  ABSL_DIE_IF_NULL(pending_run->done)->WaitForNotification();

  // Remove our progress callback.
  if (my_progress_callback.has_value()) {
    absl::MutexLock l(&queue_lock_);
    progress_callbacks_.erase(my_progress_callback.value());
  }

  return pending_run.status;
}

template<typename Key>
void Reanimator<Key>::Loop() {
  auto queue_not_empty = [this] {
    queue_lock_.AssertReaderHeld();
    return !queue_.empty();
  };

  absl::Mutex l(&queue_lock_);
  for (;;) {
    lock_.Await(absl::Condition(&queue_not_empty));
    // Run the action at the end of the queue (the one with the greatest key).
    // Note that the run is still in the queue when the action is executed.
    auto const rbegin = queue_.rbegin();
    auto const& [key, handle] = *queue_.rbegin();

    queue_lock_.Unlock();
    handle->status = action_(key);
    queue_lock_.Lock();

    // Run the progress callbacks.  This happens before unblocking the waiters
    // so that the waiters see their own progress.
    for (auto const& progress_callback : progress_callbacks_) {
      progress_callback(key, handle->status);
    }

    // Unblock the waiters.
    if (handle->done != nullptr) {
      handle->done->Notify();
    }

    // Remove the pending run from the queue.  It will be destroyed once all its
    // handles are gone.
    queue_.erase(rbegin);
  }
}

}  // namespace internal
}  // namespace _reanimator
}  // namespace base
}  // namespace principia
