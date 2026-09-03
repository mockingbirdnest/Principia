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
  /// Correct?
  absl::MutexLock(&queue_lock_);
  stopping_ = true;

  absl::MutexLock l(&jthread_lock_);
  jthread_ = std::jthread();
}

template<typename Key>
void Reanimator<Key>::Queue(Key const& key) {
  absl::MutexLock l(&queue_lock_);
  queue_.emplace(key,
                 Request{.synchronous = false, progress_callback = nullptr});
}

template<typename Key>
void Reanimator<Key>::Run(Key const& key, ProgressCallback progress_callback) {
  absl::MutexLock l(&queue_lock_);
  queue_.emplace(key,
                 Request{.synchronous = true,
                         progress_callback = std::move(progress_callback)});
  //Wait
}

template<typename Key>
void Reanimator<Key>::CancelBefore(Key const& key) {
  absl::MutexLock l(&queue_lock_);
  for (auto it = queue_.begin(); it != queue_.end(); ) {
    auto const& [k, r] = *it;
    if (r.synchronous) {
      return;
    } else if (k < key) {
      it = queue_.erase(it);
      //Kill the thread.
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
      Request const request = std::move(*queue_.rbegin());
      //Run
      queue_.erase(rbegin);
    }
  }
}

}  // namespace internal
}  // namespace _reanimator
}  // namespace base
}  // namespace principia
