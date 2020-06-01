
#pragma once

#include "numerics/pid.hpp"

#include <numeric>

#include "numerics/finite_difference.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_pid {

namespace si = quantities::si;

template<typename Input,
         typename Output,
         int horizon,
         int finite_difference_order>
PID<Input, Output, horizon, finite_difference_order>::PID(Kp const& kp,
                                                          Ki const& ki,
                                                          Kd const& kd)
    : kp_(kp), ki_(ki), kd_(kd) {}

template<typename Input,
         typename Output,
         int horizon,
         int finite_difference_order>
void PID<Input, Output, horizon, finite_difference_order>::Clear() {
  errors_.clear();
}

template<typename Input,
         typename Output,
         int horizon,
         int finite_difference_order>
Output PID<Input, Output, horizon, finite_difference_order>::
ComputeControlVariable(Input const& process_variable,
                       Input const& set_point,
                       Instant const& t) {
  // The PID is not prepared to receive inputs that are not equally spaced.
  // Clear it if that happens, but be very tolerant because KSP uses single-
  // precision float.
  Time Δt;
  if (previous_t_.has_value()) {
    Δt = t - previous_t_.value();
    if (previous_Δt_.has_value() && (Δt < previous_Δt_.value() * (1 - 1e-3) ||
                                     Δt > previous_Δt_.value() * (1 + 1e-3))) {
      LOG(WARNING) << "Resetting PID, " << Δt << " != " << previous_Δt_.value();
      Clear();
    }
    previous_Δt_ = Δt;
  }
  previous_t_ = t;

  errors_.push_back(set_point - process_variable);
  int size = errors_.size();
  if (size <= horizon) {
    // If we don't have enough error history, just return the latest error.
    return si::Unit<Kp> * errors_.back();
  }
  errors_.pop_front();
  --size;

  CHECK_LT(Time(), Δt);

  // The last error.
  Difference<Input> const proportional = errors_.back();

  // Compute the integral of the errors.  This is the dumbest possible algorithm
  // because we do not care too much about the accuracy.
  Product<Difference<Input>, Time> const integral =
      std::accumulate(errors_.cbegin(), errors_.cend(), Difference<Input>{}) *
      Δt;

  // Compute the derivative of the errors.
  std::array<Difference<Input>, finite_difference_order> errors{};
  auto it = errors_.cend() - finite_difference_order;
  for (int i = 0; i < finite_difference_order; ++i, ++it) {
    errors[i] = *it;
  }
  Variation<Input> const derivative =
      FiniteDifference(errors, Δt, finite_difference_order - 1);

  return kp_ * proportional + ki_ * integral + kd_ * derivative;
}

}  // namespace internal_pid
}  // namespace numerics
}  // namespace principia
