#pragma once

#include <functional>
#include <thread>

#include "absl/synchronization/mutex.h"

namespace principia {
namespace base {

// A minimal implementation of the C++20 jthread library, intended to be
// compatible.

class stop_source;

// Provides the means to check if a stop request has been made or can be made,
// for its associated stop_source object.
// https://en.cppreference.com/w/cpp/thread/stop_token
class stop_token {
 public:
  bool stop_requested() const;

 private:
  explicit stop_token(stop_source* ss);

  stop_source& get_stop_source() const;

  stop_source* const ss_;

  friend class stop_callback;
  friend class stop_source;
};

// Provides the means to issue a stop request. A stop request made for one
// stop_source object is visible to all stop_sources and stop_tokens of the same
// associated stop-state.
// https://en.cppreference.com/w/cpp/thread/stop_source
class stop_source {
 public:
  stop_source();

  bool request_stop();

  bool stop_requested() const;

  stop_token get_token() const;

 private:
  using Callback = std::function<void()>;

  void Register(Callback callback);

  absl::Mutex lock_;
  bool stopped_ = false;
  std::vector<Callback> callbacks_;

  friend class stop_callback;
};

// An RAII object type that registers a callback function for an associated
// stop_token object, such that the callback function will be invoked when the
// stop_token's associated stop_source is requested to stop.
// https://en.cppreference.com/w/cpp/thread/stop_callback
class stop_callback {
 public:
  explicit stop_callback(stop_token const& st, std::function<void()> callback);
  ~stop_callback();
};

// A single thread of execution which and can be cancelled/stopped in certain
// situations.
// https://en.cppreference.com/w/cpp/thread/jthread
class jthread {
 public:
  jthread();

  template<typename... Args>
  explicit jthread(std::function<void(stop_token const&, Args...)>,
                   Args&&... args);

  void join();

  void detach();

  bool request_stop();

  stop_source get_stop_source() const;

  stop_token get_stop_token() const;

 private:
  std::thread thread_;
  stop_source stop_source_;
};

}  // namespace base
}  // namespace principia