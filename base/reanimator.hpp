#pragma once

#include <functional>
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

template<typename Key>
class Reanimator {
 public:
  using Task = std::function<absl::Status(Key const&)>;
  using ProgressCallback =
      std::function<void(Key const& key, absl::Status const& status)>;

  explicit Reanimator(Task task);

  //Start is idempotent.
  void Start();
  void Stop();

  //Asynchronous.
  void Queue(Key const& key);

  //Synchronous.
  void Run(Key const& key);
  void Run(Key const& key, ProgressCallback progress_callback);

  //Cancels < before_key.
  void Cancel(Key const& before_key);  //Beware of locking.

 private:
  struct Request {
    absl::Notification* done = nullptr;
    ProgressCallback progress_callback = nullptr;
  };

  void Loop();

  Task const task_;

  absl::Mutex queue_lock_;
  bool stopping_ = false ABSL_GUARDED_BY(queue_lock_);
  //Increasing key.
  absl::btree_multimap<Key, Request> queue_ ABSL_GUARDED_BY(queue_lock_);

  absl::Mutex jthread_lock_;
  std::jthread jthread_ ABSL_GUARDED_BY(queue_lock_);
};

}  // namespace internal

using internal::Reanimator;

}  // namespace _reanimator
}  // namespace base
}  // namespace principia

#include "base/reanimator_body.hpp"