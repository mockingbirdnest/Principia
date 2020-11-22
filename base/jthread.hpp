#pragma once

#include <functional>
#include <memory>
#include <set>
#include <thread>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"

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
  bool stop_requested() const;

 private:
  explicit stop_token(not_null<StopState*> stop_state);

  StopState& get_stop_state() const;

  not_null<StopState*> const stop_state_;

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

  friend class StopState;
};

// A single thread of execution which and can be cancelled/stopped in certain
// situations.
// https://en.cppreference.com/w/cpp/thread/jthread
class jthread {
 public:
  jthread() = default;

  template<typename... Args>
  explicit jthread(std::function<void(stop_token const&, Args...)> f,
                   Args&&... args);

  void join();

  void detach();

  bool request_stop();

  stop_source get_stop_source() const;

  stop_token get_stop_token() const;

 private:
  std::unique_ptr<StopState> stop_state_;
  std::thread thread_;
};

}  // namespace internal_jthread

using internal_jthread::jthread;
using internal_jthread::stop_callback;
using internal_jthread::stop_source;
using internal_jthread::stop_token;

}  // namespace base
}  // namespace principia
