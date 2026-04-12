#pragma once

#include "base/stoppable_thread.hpp"

#include <memory>
#include <set>
#include <stop_token>
#include <thread>
#include <utility>

namespace principia {
namespace base {
namespace _stoppable_thread {
namespace internal {

template<typename Function, typename... Args>
std::jthread MakeStoppableThread(Function&& f, Args&&... args) {
  return std::jthread(
      [f](std::stop_token const& st, Args&&... args) {
        // This assignment happens on the thread of the jthread.
        this_stoppable_thread::stop_token_ = st;
        f(std::forward<Args>(args)...);
      },
      std::forward<Args>(args)...);
}

inline std::stop_token this_stoppable_thread::get_stop_token() {
  return stop_token_;
}

}  // namespace internal
}  // namespace _stoppable_thread
}  // namespace base
}  // namespace principia
