#pragma once

#include <chrono>
#include <functional>
#include <optional>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/synchronization/mutex.h"
#include "base/jthread.hpp"
#include "base/macros.hpp"

namespace principia {
namespace base {
namespace internal_recurring_thread {

// A stoppable thread that supports cyclical execution of an action.  It is
// connected to two monodirectional channels that can (optionally) hold a value
// of |Input| (for incoming data) or |Output| (for outgoing data), respectively.
// The action is run to transform the input into the output.  This class and its
// subclasses are thread-safe.  The base class is used to factor code common to
// the various template specializations and should not be used directly.
class BaseRecurringThread {
 public:
  virtual ~BaseRecurringThread() = default;

  // Starts or stops the thread.  These functions are idempotent.  Note that the
  // thread is also stopped by the destruction of this object.
  void Start();
  void Stop();

  // Stop followed by Start, atomically.
  void Restart();

 protected:
  // Constructs a stoppable thread that runs no more frequently than at the
  // specified |period| (and less frequently if no input was provided).  At
  // construction the thread is in the stopped state.
  explicit BaseRecurringThread(std::chrono::milliseconds period);

  // Repeatedly calls RunAction no more frequently than at the specified period.
  absl::Status RepeatedlyRunAction();

  // Overidden by subclasses to actually run the action.
  virtual absl::Status RunAction() = 0;

 private:
  std::chrono::milliseconds const period_;

  absl::Mutex jthread_lock_;
  jthread jthread_ GUARDED_BY(jthread_lock_);
};

// A template for an action that returns a value.
template<typename Input, typename Output = void>
class RecurringThread : public BaseRecurringThread {
 public:
  // If an action returns an error, no output in written to the output channel.
  using Action = std::function<absl::StatusOr<Output>(Input)>;

  // Constructs a stoppable thread that executes the given |action| no more
  // frequently than at the specified |period| (and less frequently if no input
  // was provided).  At construction the thread is in the stopped state.
  RecurringThread(Action action,
                  std::chrono::milliseconds period);

  // Overwrites the contents of the input channel.  The |input| data will be
  // either picked by the next execution of |action|, or overwritten by the next
  // call to |Put|.
  void Put(Input input);

  // Extracts data from the output channel, if there is any.
  std::optional<Output> Get();

 private:
  absl::Status RunAction() override;

  Action const action_;

  absl::Mutex input_output_lock_;
  std::optional<Input> input_ GUARDED_BY(input_output_lock_);
  std::optional<Output> output_ GUARDED_BY(input_output_lock_);
};

// A template for an action that returns no value.
template<typename Input>
class RecurringThread<Input, void> : public BaseRecurringThread {
 public:
  using Action = std::function<absl::Status(Input)>;

  // Constructs a stoppable thread that executes the given |action| no more
  // frequently than at the specified |period| (and less frequently if no input
  // was provided).  At construction the thread is in the stopped state.
  RecurringThread(Action action,
                  std::chrono::milliseconds period);

  // Overwrites the contents of the input channel.  The |input| data will be
  // either picked by the next execution of |action|, or overwritten by the next
  // call to |Put|.
  void Put(Input input);

 private:
  absl::Status RunAction() override;

  Action const action_;

  absl::Mutex input_lock_;
  std::optional<Input> input_ GUARDED_BY(input_lock_);
};

}  // namespace internal_recurring_thread

using internal_recurring_thread::RecurringThread;

}  // namespace base
}  // namespace principia

#include "base/recurring_thread_body.hpp"
