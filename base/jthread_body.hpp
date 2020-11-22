#pragma once

#include "base/jthread.hpp"

#include "base/macros.hpp"

namespace principia {
namespace base {
namespace internal_jthread {

class StopState {
 public:
  StopState() = default;

  bool request_stop();

  bool stop_requested() const;

  void Register(not_null<stop_callback*> callback);
  void Unregister(not_null<stop_callback*> callback);

 private:
  mutable absl::Mutex lock_;
  bool stop_requested_ GUARDED_BY(lock_) = false;
  std::set<not_null<stop_callback*>> callbacks_ GUARDED_BY(lock_);
};

inline bool StopState::request_stop() {
  // NOTE(phl): If performance matters here we could do double-locking.
  std::set<not_null<stop_callback*>> callbacks;
  {
    absl::MutexLock l(&lock_);
    if (stop_requested_) {
      return false;
    } else {
      stop_requested_ = true;
      callbacks.swap(callbacks_);
    }
  }
  for (auto const callback : callbacks_) {
    callback->Run();
  }
  return true;
}

inline bool StopState::stop_requested() const {
  absl::ReaderMutexLock l(&lock_);
  return stop_requested_;
}

inline void StopState::Register(not_null<stop_callback*> const callback) {
  absl::MutexLock l(&lock_);
  callbacks_.insert(callback);
}

inline void StopState::Unregister(not_null<stop_callback*> const callback) {
  absl::MutexLock l(&lock_);
  callbacks_.erase(callback);
}

inline bool stop_token::stop_requested() const {
  return stop_state_->stop_requested();
}

inline stop_token::stop_token(not_null<StopState*> const stop_state)
    : stop_state_(stop_state) {}

inline StopState& stop_token::get_stop_state() const {
  return *stop_state_;
}

inline bool stop_source::request_stop() {
  return stop_state_->request_stop();
}

inline bool stop_source::stop_requested() const {
  return stop_state_->stop_requested();
}

inline stop_token stop_source::get_token() const {
  return stop_token(stop_state_);
}

inline stop_source::stop_source(not_null<StopState*> const stop_state)
    : stop_state_(stop_state) {}

inline stop_callback::stop_callback(stop_token const& st,
                             std::function<void()> callback)
    : callback_(std::move(callback)) {
  st.get_stop_state().Register(this);
}

inline stop_callback::~stop_callback() {}

inline void stop_callback::Run() const {
  callback_();
}

template<typename Function, typename... Args>
jthread::jthread(Function&& f, Args&&... args)
    : stop_state_(std::make_unique<StopState>()),
      thread_(std::move(f),
              stop_token(stop_state_.get()),
              std::forward<Args>(args)...) {}

inline jthread::~jthread() {
  stop_state_->request_stop();
  if (thread_.joinable()) {
    thread_.join();
  }
}

inline void jthread::join() {
  thread_.join();
}

inline void jthread::detach() {
  thread_.detach();
}

inline bool jthread::request_stop() {
  return stop_state_->request_stop();
}

inline stop_source jthread::get_stop_source() const {
  return stop_source(stop_state_.get());
}

inline stop_token jthread::get_stop_token() const {
  return stop_token(stop_state_.get());
}

}  // namespace internal_jthread
}  // namespace base
}  // namespace principia