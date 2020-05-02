
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

  // Adds the error between the two values to the state of the PID and
  // returns an apparent value derived from the actual value and the control
  // variable.
  Value ComputeValue(Value const& apparent,
                     Value const& actual,
                     Time const& Δt);

private:
  double const kp_;
  Inverse<Time> const ki_;
  Time const kd_;

  // The front element is the oldest.
  std::deque<Value> errors_;
};

}  // namespace internal_pid

using internal_pid::PID;

}  // namespace numerics
}  // namespace principia

#include "numerics/pid_body.hpp"
