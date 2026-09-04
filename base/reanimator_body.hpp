#pragma once

#include "base/reanimator.hpp"

#include "absl/log/check.h"
#include "absl/log/log.h"
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

  auto queue_empty = [this] {
    queue_lock_.AssertReaderHeld();
    LOG(ERROR) << "Stop condition E: " << queue_.empty();
    return queue_.empty();
  };

  // Mark the thread as stopped and wait for all the calls to complete.  This
  // includes the best-effort calls.  The actions that have a
  // `RETURN_IF_STOPPED` should finish quickly.  Note that after this point it's
  // not possible to enqueue new runs.
  {
    absl::MutexLock l(&queue_lock_);
    LOG(ERROR) << "Stop request_stop";
    jthread_.request_stop();
    stopped_ = true;
    queue_lock_.Await(absl::Condition(&queue_empty));
  }

  LOG(ERROR) << "Stop joining";
  jthread_.join();
  LOG(ERROR) << "Stop joined";
}

template<typename Key>
void Reanimator<Key>::RunBestEffort(Key const& key) {
  absl::MutexLock l(&queue_lock_);
  CHECK(!stopped_);
  queue_.emplace(key, std::make_shared<PendingRun>());
}

template<typename Key>
typename Reanimator<Key>::Handle Reanimator<Key>::RunGuaranteed(
    Key const& key) {
  absl::MutexLock l(&queue_lock_);
  CHECK(!stopped_);

  // Queue the call.  The `Notification` will be notified once it has run.
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
    absl::MutexLock l(&queue_lock_);
    for (auto it = queue_.begin(); it != queue_.end();) {
      auto const& [key, run] = *it;
      if (run->done != nullptr) {
        return;
      } else if (key < before_key) {
        it = queue_.erase(it);
        LOG(ERROR) << "Cancel queue_size: " << queue_.size();
        // If the queue is now empty, we need to mark the thread as stopped.  If
        // there is a running action with a `RETURN_IF_STOPPED`, it should
        // complete quickly.
        if (it == queue_.end()) {
          LOG(ERROR) << "Cancel request_stop";
          jthread_.request_stop();
        }
      } else {
        return;
      }
    }
  }

  LOG(ERROR) << "Cancel recreating thread";
  jthread_.join();
  jthread_ = MakeStoppableThread([this]() { RepeatedRunActions(); });
  LOG(ERROR) << "Cancel recreated thread";
}

template<typename Key>
absl::Status Reanimator<Key>::Wait(Handle const handle,
                                   ProgressCallback progress_callback) {
  // This object won't go away since we hold `handle`.
  auto const& pending_run = *ABSL_DIE_IF_NULL(handle);

  // Insert my progress callback, if any, and retain an iterator on it.
  std::optional<ProgressCallbacks::const_iterator> my_progress_callback;
  if (progress_callback != nullptr) {
    absl::MutexLock l(&queue_lock_);
    progress_callbacks_.push_front(std::move(progress_callback));
    my_progress_callback = progress_callbacks_.begin();
  }

  ABSL_DIE_IF_NULL(pending_run.done)->WaitForNotification();

  // Remove my progress callback.
  if (my_progress_callback.has_value()) {
    absl::MutexLock l(&queue_lock_);
    progress_callbacks_.erase(my_progress_callback.value());
  }

  return pending_run.status;
}

template<typename Key>
void Reanimator<Key>::RepeatedRunActions() {
  auto queue_not_empty_or_stop_requested = [this] {
    queue_lock_.AssertReaderHeld();
    // All requests to stop must happen under `queue_lock_`.  Don't use
    // `this_stoppable_thread` here, the condition can be evaluated on any
    // thread.
    LOG(ERROR) << "RepeatedRunActions condition E: " << queue_.empty()
               << " SR: " << jthread_.get_stop_token().stop_requested();
    return !queue_.empty() || jthread_.get_stop_token().stop_requested();
  };

  absl::MutexLock l(&queue_lock_);
  for (;;) {
    queue_lock_.Await(absl::Condition(&queue_not_empty_or_stop_requested));

    LOG(ERROR) << "RepeatedRunActions awaited E: " << queue_.empty() << " SR: "
               << this_stoppable_thread::get_stop_token().stop_requested();
    if (queue_.empty()) {
      CHECK(this_stoppable_thread::get_stop_token().stop_requested());
      LOG(ERROR) << "RepeatedRunActions breaking";
      break;
    }

    // Run the action at the end of the queue (the one with the greatest key).
    // Note that the run is still in the queue when the action is executed.  It
    // is important to retain a copy of `handle` because the queue might change
    // (due to `Cancel`) while we run the action.
    auto const [key, handle] = *queue_.rbegin();

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
