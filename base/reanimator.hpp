#pragma once

#include <functional>
#include <thread>

#include "absl/container/btree_map.h"
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

  void Start();
  void Stop();

  //Asynchronous.
  void Queue(Key const& key);

  //Synchronous.
  void Run(Key const& key, ProgressCallback progress_callback);

  void CancelBefore(Key const& key);  //Beware of locking.

 private:
  struct Request {
    bool synchronous = true;
  };

  Task const task_;
  std::jthread thread_;
  absl::btree_multimap<Key, Request> queue_;
};

}  // namespace internal

}  // namespace _reanimator
}  // namespace base
}  // namespace principia
