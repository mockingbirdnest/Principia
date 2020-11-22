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

  void Register(stop_callback* callback);
  void Unregister(stop_callback* callback);

 private:
  mutable absl::Mutex lock_;
  bool stop_requested_ GUARDED_BY(lock_) = false;
  std::set<stop_callback*> callbacks_ GUARDED_BY(lock_);
};

bool StopState::request_stop() {
  // NOTE(phl): If performance matters here we could do double-locking.
  std::set<stop_callback*> callbacks;
  {
    absl::MutexLock l(&lock_);
    if (stop_requested_) {
      return false;
    } else {
      stop_requested_ = true;
      callbacks.swap(callbacks_);
    }
  }
  for (auto* const callback : callbacks_) {
    callback->Run();
  }
  return true;
}

bool StopState::stop_requested() const {
  absl::ReaderMutexLock l(&lock_);
  return stop_requested_;
}

void StopState::Register(stop_callback* const callback) {
  absl::MutexLock l(&lock_);
  callbacks_.insert(callback);
}

void StopState::Unregister(stop_callback* const callback) {
  absl::MutexLock l(&lock_);
  callbacks_.erase(callback);
}

bool stop_token::stop_requested() const {
  return stop_state_->stop_requested();
}

stop_token::stop_token(StopState* const stop_state) : stop_state_(stop_state) {}

StopState& stop_token::get_stop_state() const {
  return *stop_state_;
}

bool stop_source::request_stop() {
  return stop_state_->request_stop();
}

bool stop_source::stop_requested() const {
  return stop_state_->stop_requested();
}

stop_token stop_source::get_token() const {
  return stop_token(stop_state_);
}

stop_callback::stop_callback(stop_token const& st,
                             std::function<void()> callback)
    : callback_(std::move(callback)) {
  st.get_stop_state().Register(this);
}

stop_callback::~stop_callback() {}

void stop_callback::Run() const {
  callback_();
}

template<typename... Args>
jthread::jthread(std::function<void(stop_token const&, Args...)> f,
                 Args&&... args)
    : thread_([f = std::move(f), args = std::move(args...)]() { f(args...); }),
      stop_state_(std::make_unique<StopState>()) {}

void jthread::join() {
  thread_.join();
}

void jthread::detach() {
  thread_.detach();
}

bool jthread::request_stop() {
  return stop_state_->request_stop();
}
stop_source jthread::get_stop_source() const {
  return stop_source(stop_state_.get());
}

stop_token jthread::get_stop_token() const {
  return stop_token(stop_state_.get());
}

}  // namespace internal_jthread
}  // namespace base
}  // namespace principia