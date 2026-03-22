#pragma once

#include <functional>
#include <stop_token>
#include <thread>

#include "absl/status/status.h"

namespace principia {
namespace base {
namespace _stoppable_thread {
namespace internal {

// `f` *does not* take a stop_token as its first parameter.
template<typename Function, typename... Args>
static std::jthread MakeStoppableThread(Function&& f, Args&&... args);

class this_stoppable_thread {
 public:
  static std::stop_token get_stop_token();

 private:
  inline static thread_local std::stop_token stop_token_;

  template<typename Function, typename... Args>
  friend std::jthread MakeStoppableThread(Function&& f, Args&&... args);
};

#define RETURN_IF_STOPPED                                             \
  do {                                                                \
    if (::principia::base::_stoppable_thread::this_stoppable_thread:: \
            get_stop_token()                                          \
                .stop_requested()) {                                  \
      return ::absl::CancelledError("Cancelled by stop token");       \
    }                                                                 \
  } while (false)

}  // namespace internal

using internal::MakeStoppableThread;
using internal::this_stoppable_thread;

}  // namespace _stoppable_thread
}  // namespace base
}  // namespace principia

#include "base/stoppable_thread_body.hpp"
