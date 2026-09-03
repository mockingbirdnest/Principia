#pragma once

#include "base/reanimator.hpp"
#include "base/stoppable_thread.hpp"

namespace principia {
namespace base {
namespace _reanimator {
namespace internal {

template<typename Key>
Reanimator<Key>::Reanimator(Task task) : task_(std::move(task)) {}

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
void Reanimator<Key>::Queue(Key const& key) {
  absl::MutexLock l(&queue_lock_);
  queue_.emplace(key, Request{});
}

template<typename Key>
void Reanimator<Key>::Run(Key const& key) {
  Run(key, /*progress_callback=*/nullptr);
}

template<typename Key>
void Reanimator<Key>::Run(Key const& key, ProgressCallback progress_callback) {
  auto const done = std::make_unique<Notification>();
  absl::MutexLock l(&queue_lock_);

  // Queue the call.  The `Notification` will be notified once it has run.
  auto const [it, _] = queue_.emplace(
      key,
      Request{.done = done.get(),
              .progress_callback = std::move(progress_callback)});

  // Wait for the queued run to execute.  Since we never cancel synchronous
  // calls we are sure that the `Notification` will ultimately be notified.
  done->WaitForNotification();
}

template<typename Key>
void Reanimator<Key>::Cancel(Key const& before_key) {
  absl::MutexLock l(&queue_lock_);
  for (auto it = queue_.begin(); it != queue_.end();) {
    auto const& [key, request] = *it;
    if (request.synchronous) {
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
      auto const [key, request] = *queue_.rbegin();
      absl::Status const status = task_(key);
      if (request.done != nullptr) {
        request.done->Notify();
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
