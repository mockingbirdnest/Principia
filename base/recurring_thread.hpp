#pragma once

#include <chrono>
#include <functional>
#include <optional>

#include "absl/synchronization/mutex.h"
#include "base/jthread.hpp"
#include "base/macros.hpp"
#include "base/status.hpp"
#include "base/status_or.hpp"

namespace principia {
namespace base {
namespace internal_recurring_thread {

// A stoppable thread that supports cyclical execution of an action.  It is
// connected to two monodirectional channels that can (optionally) hold a value
// of |Input| (for incoming data) or |Output| (for outgoing data), respectively.
// The action is run to transform the input into the output.  This class is
// thread-safe.
template<typename Input, typename Output>
class RecurringThread {
 public:
  using Action = std::function<StatusOr<Output>(Input)>;

  // Constructs a stoppable thread that executes the given |action| no more
  // frequently than at the specified |period| (and less frequently if no input
  // was provided).  At construction the thread is in the stopped state.
  RecurringThread(Action action,
                  std::chrono::milliseconds period);

  // Starts or stops the thread.  These functions are idempotent.  Note that the
  // thread is also stopped by the destruction of this object.
  void Start();
  void Stop();

  // Overwrites the contents of the input channel.  The |input| data will be
  // either picked by the next execution of |action|, or overwritten by the next
  // call to |Put|.
  void Put(Input input);

  // Extracts data from the output channel, if there is any.
  std::optional<Output> Get();

 private:
  Status RepeatedlyRunAction();

  Action const action_;
  std::chrono::milliseconds const period_;

  absl::Mutex jthread_lock_;
  jthread jthread_ GUARDED_BY(jthread_lock_);

  absl::Mutex input_output_lock_;
  std::optional<Input> input_ GUARDED_BY(input_output_lock_);
  std::optional<Output> output_ GUARDED_BY(input_output_lock_);
};

}  // namespace internal_recurring_thread

using internal_recurring_thread::RecurringThread;

}  // namespace base
}  // namespace principia

#include "base/recurring_thread_body.hpp"
