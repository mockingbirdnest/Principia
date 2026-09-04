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
  // Should wait for the sync calls.
  absl::MutexLock(&queue_lock_);
  stopping_ = true;

  absl::MutexLock l(&jthread_lock_);
  jthread_ = std::jthread();
}

template<typename Key>
void Reanimator<Key>::RunAsynchronously(Key const& key) {
  absl::MutexLock l(&queue_lock_);
  queue_.emplace(key, PendingRun{});
}

template<typename Key>
typename Reanimator<Key>::Handle Reanimator<Key>::RunSynchronously(
    Key const& key) {
  absl::MutexLock l(&queue_lock_);

  // Queue the call.  The `Notification` will be notified once it has run.
  auto const [it, _] = queue_.emplace(
      key,
      PendingRun{.done = std::make_unique<Notification>()});

  return &it->second;
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
      // Kill the thread.
    } else {
      return;
    }
  }
}

template<typename Key>
absl::Status Reanimator<Key>::Wait(Handle const handle,
                                   ProgressCallback progress_callback) {
  auto const& pending_run = *ABSL_DIE_IF_NULL(handle);
  ABSL_DIE_IF_NULL(pending_run->done)->WaitForNotification();
  //Progress, status
}

template<typename Key>
void Reanimator<Key>::Loop() {
  auto queue_not_empty_or_stopping = [this] {
    queue_lock_.AssertReaderHeld();
    return !queue_.empty() || stopping_;
  };

  absl::Mutex l(&queue_lock_);
  for (;;) {
    lock_.Await(absl::Condition(&queue_not_empty_or_stopping));
    if (stopping_) {
      return;
    } else {
      auto const rbegin = queue_.rbegin();
      auto const [key, run] = *queue_.rbegin();
      absl::Status const status = action_(key);
      if (run.done != nullptr) {
        run.done->Notify();
      }
      // Progress
      queue_.erase(rbegin);
    }
  }
}

}  // namespace internal
}  // namespace _reanimator
}  // namespace base
}  // namespace principia
