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

template<typename Key>
class Reanimator {
  struct PendingRun;
  using Handle = std::shared_ptr<PendingRun>;

 public:
  using Action = std::function<absl::Status(Key const&)>;
  using ProgressCallback =
      std::function<void(Key const& key, absl::Status const& status)>;

  explicit Reanimator(Action action);

  ///Start is idempotent.
  void Start();
  void Stop();

  void RunAsynchronously(Key const& key);

  Handle RunSynchronously(Key const& key);

  ///Cancels < before_key.
  void Cancel(Key const& before_key);

  absl::Status Wait(Handle handle,
                    ProgressCallback progress_callback = nullptr);

 private:
  ///Iterator stability.
  using ProgressCallbacks = std::list<ProgressCallback>;

  struct PendingRun {
    std::unique_ptr<absl::Notification> done;//Synch
    absl::Status status;
  };

  void Loop();

  Action const action_;

  absl::Mutex queue_lock_;
  bool stopping_ = false ABSL_GUARDED_BY(queue_lock_);
  ProgressCallbacks progress_callbacks_ ABSL_GUARDED_BY(queue_lock_);
  ///Increasing key.  Pointer stability.
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