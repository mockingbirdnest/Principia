#pragma once

#include <deque>
#include <functional>
#include <future>
#include <list>
#include <thread>

#include "absl/base/thread_annotations.h"
#include "absl/synchronization/mutex.h"

namespace principia {
namespace base {

// A pool of threads that are created at construction and to which functions can
// be added for asynchronous execution.  This class is thread-safe.
template<typename T>
class ThreadPool final {
 public:
  // Constructs a pool with the given number of threads.
  explicit ThreadPool(std::int64_t pool_size);

  ~ThreadPool();

  // Adds a call to the execution queue, and returns a future that the client
  // may use to wait until execution of |function| has completed and to extract
  // the result.
  std::future<T> Add(std::function<T()> function);

 private:
  // The queue element contains a |function| to execute and a |promise| used to
  // communicate the result to the caller.
  struct Call {
    std::function<T()> function;
    std::promise<T> promise;
  };

  // The loop executed on each thread to extract an element from the queue,
  // execute it, and set its result in the promise.
  void DequeueCallAndExecute();

  absl::Mutex lock_;
  bool shutdown_ GUARDED_BY(lock_) = false;
  std::deque<Call> calls_ GUARDED_BY(lock_);

  std::list<std::thread> threads_;
};

}  // namespace base
}  // namespace principia

#include "base/thread_pool_body.hpp"
