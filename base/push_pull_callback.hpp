#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <thread>
#include <tuple>

#include "absl/synchronization/mutex.h"
#include "base/macros.hpp"  // 🧙 For GUARDED_BY.

namespace principia {
namespace base {
namespace _push_pull_callback {
namespace internal {

// A helper class for callbacks from unmanaged code (C++) to managed code (C#).
// While we could support true callbacks in the generator, it is difficult to
// implement journalling because a journal entry would have to include all the
// parameters passed to callbacks and the returned results.
// Instead we invert the flow of control (just like we do for serialization) and
// expect the managed code to pull the arguments and push the result of the
// computations that it is executing.
// Note that calling the managed and unmanaged APIs from the same thread will
// inevitably cause deadlocks.  See `PushPullExecutor` below for a solution to
// this.
template<typename Result, typename... Arguments>
class PushPullCallback {
 public:
  // The managed API, called to extract the arguments for the unmanaged callback
  // and return its result.  `Pull` returns false if there are no more arguments
  // to be processed and the managed code should stop its iteration.
  bool Pull(Arguments&... arguments);
  void Push(Result result);

  // Used on the unmanaged side to use this object in a context that requires a
  // function.
  std::function<Result(Arguments...)> ToStdFunction();

  // Used on the unmanaged side to indicate that the computation has finished.
  // After a call to this method, `Pull` always returns false.
  void Shutdown();

 private:
  // The unmanaged API, called by the function returned by `ToStdFunction`.
  void Push(Arguments... arguments);
  Result Pull();

  // These functions return a (held) `MutexLock` that the caller should use to
  // ensure proper release of `lock_`.
  std::unique_ptr<absl::MutexLock> WaitUntilHasArgumentsOrShuttingDownAndLock();
  std::unique_ptr<absl::MutexLock> WaitUntilHasResultAndLock();

  absl::Mutex lock_;
  std::optional<std::tuple<Arguments...>> arguments_ GUARDED_BY(lock_);
  std::optional<Result> result_ GUARDED_BY(lock_);
  bool shutdown_ GUARDED_BY(lock_) = false;
};

// A helper class to execute a task that takes a callback from unmanaged code to
// managed code and returns a value of type `T`.  The task is executed on a
// separate thread, so calls to the`PushPullCallback` don't cause deadlocks.
template<typename T,
         typename Result, typename... Arguments>
class PushPullExecutor {
 public:
  using Task = std::function<T(std::function<Result(Arguments...)>)>;

  explicit PushPullExecutor(Task task);
  ~PushPullExecutor();

  // Returns the internal `PushPullCallback` object that is used by the managed
  // code to pull arguments and push results.
  PushPullCallback<Result, Arguments...>& callback();
  PushPullCallback<Result, Arguments...> const& callback() const;

  // Returns the result of the task passed at construction.  Must only be called
  // once `PushPullCallback::Pull` has indicated that the task has finished.
  T get();

 private:
  PushPullCallback<Result, Arguments...> callback_;
  mutable absl::Mutex lock_;
  std::optional<T> result_ GUARDED_BY(lock_);

  // This must come last as it references the other member variables, see #4136.
  std::thread thread_;
};

}  // namespace internal

using internal::PushPullCallback;
using internal::PushPullExecutor;

}  // namespace _push_pull_callback
}  // namespace base
}  // namespace principia

#include "base/push_pull_callback_body.hpp"
