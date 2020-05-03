
#pragma once

#include "numerics/pid.hpp"

#include <numeric>

#include "numerics/finite_difference.hpp"
#include "numerics/ulp_distance.hpp"

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
    Value const& process_variable,
    Value const& set_point,
    Time const& Δt) {
  // The PID is not prepared to receive inputs that are not equally spaced.
  // Clear it if that happens, but be very tolerant because KSP uses single-
  // precision float.
  if (previous_Δt_.has_value() && (Δt < previous_Δt_.value * (1 - 1e-3) ||
                                   Δt > previous_Δt_.value * (1 + 1e-3))) {
    Clear();
  }
  previous_Δt_ = Δt;

  errors_.push_back(process_variable - set_point);
  int size = errors_.size();
  if (size <= horizon) {
    // If we don't have enough error history, just return our input.
    return process_variable;
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

  return set_point + (kp_ * proportional + ki_ * integral + kd_ * derivative);
}

}  // namespace internal_pid
}  // namespace numerics
}  // namespace principia
