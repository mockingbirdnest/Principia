
#pragma once

#include <deque>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_pid {

using quantities::Inverse;
using quantities::Time;

template<typename Value, int horizon, int finite_difference_order>
class PID {
public:
  static_assert(finite_difference_order <= horizon);

  PID(double kp, Inverse<Time> const& ki, Time const& kd);

  // Clears the state of the PID.
  void Clear();

  // Adds the error between the process variable and the set-point to the state
  // of the PID and returns an updated process variable derived from the
  // set-point and the control variable.
  Value ComputeValue(Value const& process_variable,
                     Value const& set_point,
                     Time const& Δt);

private:
  double const kp_;
  Inverse<Time> const ki_;
  Time const kd_;

  // The PID is not prepared to receive inputs that are not equally spaced.
  std::optional<Time> previous_Δt_;

  // The front element is the oldest.
  std::deque<Value> errors_;
};

}  // namespace internal_pid

using internal_pid::PID;

}  // namespace numerics
}  // namespace principia

#include "numerics/pid_body.hpp"
