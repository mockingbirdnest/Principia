#pragma once

#include <atomic>
#include <functional>
#include <queue>
#include <mutex>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

class Bundle {
  Bundle(int threads);

  using Task = std::function<Status()>;
  void Add(Task task);
  Status Join();

 private:
  std::mutex lock_;
  std::condition_variable tasks_not_empty_or_joining_;

  bool joining_ GUARDED_BY(lock_);
  std::queue<Task> tasks_ GUARDED_BY(lock_);
  Status status_ GUARDED_BY(lock_);

  std::vector<std::thread> threads_;
};

}  // namespace base
}  // namespace principia
