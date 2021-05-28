#pragma once

#include <functional>
#include <memory>
#include <set>
#include <thread>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "absl/status/status.h"

namespace principia {
namespace base {
namespace internal_jthread {

// A minimal implementation of the C++20 jthread library, intended to be
// compatible.

class StopState;

// Provides the means to check if a stop request has been made or can be made,
// for its associated stop_source object.
// https://en.cppreference.com/w/cpp/thread/stop_token
class stop_token {
 public:
  stop_token() = default;

  bool stop_requested() const;

 private:
  explicit stop_token(not_null<StopState*> stop_state);

  StopState& get_stop_state() const;

  StopState* stop_state_ = nullptr;

  friend class jthread;
  friend class stop_callback;
  friend class stop_source;
};

// Provides the means to issue a stop request. A stop request made for one
// stop_source object is visible to all stop_sources and stop_tokens of the same
// associated stop-state.
// https://en.cppreference.com/w/cpp/thread/stop_source
class stop_source {
 public:
  bool request_stop();

  bool stop_requested() const;

  stop_token get_token() const;

 private:
  explicit stop_source(not_null<StopState*> stop_state);

  not_null<StopState*> const stop_state_;

  friend class jthread;
};

// An RAII object type that registers a callback function for an associated
// stop_token object, such that the callback function will be invoked when the
// stop_token's associated stop_source is requested to stop.
// https://en.cppreference.com/w/cpp/thread/stop_callback
class stop_callback {
 public:
  explicit stop_callback(stop_token const& st, std::function<void()> callback);
  ~stop_callback();

 private:
  void Run() const;

  std::function<void()> const callback_;
  stop_token const stop_token_;

  friend class StopState;
};

// A single thread of execution which and can be cancelled/stopped in certain
// situations.
// https://en.cppreference.com/w/cpp/thread/jthread
class jthread {
 public:
  jthread() = default;

  template<typename Function, typename... Args>
  jthread(Function&& f, Args&&... args);

  jthread(jthread&& other);
  jthread& operator=(jthread&& other);

  ~jthread();

  bool joinable() const;

  void join();

  void detach();

  bool request_stop();

  stop_source get_stop_source() const;

  stop_token get_stop_token() const;

 private:
  std::unique_ptr<StopState> stop_state_;
  std::thread thread_;
};

// |f| *does not* take a stop_token as its first parameter.
template<typename Function, typename... Args>
static jthread MakeStoppableThread(Function&& f, Args&&... args);

class this_stoppable_thread {
 public:
  static stop_token get_stop_token();

 private:
  inline static thread_local stop_token stop_token_;

  template<typename Function, typename... Args>
  friend jthread MakeStoppableThread(Function&& f, Args&&... args);
};

#define RETURN_IF_STOPPED                                                   \
  do {                                                                      \
    if (::principia::base::this_stoppable_thread::get_stop_token()          \
            .stop_requested()) {                                            \
      return ::principia::base::Status(::principia::base::Error::CANCELLED, \
                                       "Cancelled by stop token");          \
    }                                                                       \
  } while (false)

}  // namespace internal_jthread

using internal_jthread::MakeStoppableThread;
using internal_jthread::jthread;
using internal_jthread::stop_callback;
using internal_jthread::stop_source;
using internal_jthread::stop_token;
using internal_jthread::this_stoppable_thread;

}  // namespace base
}  // namespace principia

#include "base/jthread_body.hpp"
