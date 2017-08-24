#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <list>
#include <mutex>
#include <queue>
#include <thread>

namespace principia {
namespace base {

template<typename T>
class ThreadPool final {
 public:
  explicit ThreadPool(std::int64_t pool_size);

  ~ThreadPool();

  std::future<T> Add(std::function<T()> function);

 private:
  struct Call {
    std::function<T()> function;
    std::promise<T> promise;
  };

  void DequeueCallAndExecute();

  std::mutex lock_;
  bool shutdown_ = false;
  std::deque<Call> calls_;
  std::condition_variable has_calls_;

  std::list<std::thread> threads_;
};

}  // namespace base
}  // namespace principia

#include "base/thread_pool_body.hpp"
