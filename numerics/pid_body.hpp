
#pragma once

#include "numerics/pid.hpp"

#include <numeric>

#include "numerics/finite_difference.hpp"

namespace principia {
namespace numerics {
namespace internal_pid {

using quantities::Product;
using quantities::Variation;

template<typename Value, int horizon, int finite_difference_order>
PID<Value, horizon, finite_difference_order>::PID(double const kp,
                                                       Inverse<Time> const& ki,
                                                       Time const& kd)
    : kp_(kp), ki_(ki), kd_(kd) {}

template<typename Value, int horizon, int finite_difference_order>
void PID<Value, horizon, finite_difference_order>::Clear() {
  errors_.clear();
}

template<typename Value, int horizon, int finite_difference_order>
Value PID<Value, horizon, finite_difference_order>::ComputeValue(
    Value const& apparent,
    Value const& actual,
    Time const& Δt) {
  errors_.push_back(apparent - actual);
  int size = errors_.size();
  if (size <= horizon) {
    // If we don't have enough error history, just return our input.
    return apparent;
  }
  errors_.pop_front();
  --size;

  // The last error.
  Value proportional = errors_.back();

  // Compute the integral of the errors.  This is the dumbest possible algorithm
  // because we do not care too much about the accuracy.
  Product<Value, Time> const integral =
      std::accumulate(errors_.cbegin(), errors_.cend(), Value{}) * Δt;

  // Compute the derivative of the errors.
  std::array<Value, finite_difference_order> errors{};
  auto it = errors_.cend() - finite_difference_order;
  for (int i = 0; i < finite_difference_order; ++i, ++it) {
    errors[i] = *it;
  }
  Variation<Value> const derivative =
      FiniteDifference(errors, Δt, finite_difference_order - 1);

  return actual + (kp_ * proportional + ki_ * integral + kd_ * derivative);
}

}  // namespace internal_pid
}  // namespace numerics
}  // namespace principia