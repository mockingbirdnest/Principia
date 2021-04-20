#pragma once

#include <chrono>
#include <functional>
#include <optional>

#include "base/jthread.hpp"

namespace principia {
namespace base {
namespace internal_recurring_thread {

template<typename Input, typename Output>
class RecurringThread {
 public:
  using Action = std::function<Output(Input)>;

  RecurringThread(Action action,
                  std::chrono::duration period);

  void Put(Input input);

  std::optional<Output> Get();

 private:
  Action const action_;
  std::chrono::duration const period_;

  jthread jthread_;

  absl::Mutex lock_;
  std::optional<Input> input_;
  std::optional<Output> output_;
};

}  // namespace internal_recurring_thread
}  // namespace base
}  // namespace principia
