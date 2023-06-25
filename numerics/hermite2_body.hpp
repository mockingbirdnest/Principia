#pragma once

#include "numerics/hermite2.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace principia {
namespace numerics {
namespace _hermite2 {
namespace internal {

template<typename Value, typename Argument>
Hermite2<Value, Argument>::Hermite2(
    std::pair<Argument, Argument> arguments,
    std::pair<Value, Value> const& values,
    Derivative1 const& derivative_first)
    : arguments_(std::move(arguments)) {
  a0_ = values.first;
  a1_ = derivative_first;
  Difference<Argument> const Δargument = arguments_.second - arguments_.first;
  // If we were given the same point twice, there is a removable singularity.
  // Otherwise, if the arguments are the same but not the values, we proceed to
  // merrily NaN away as we should.
  if (Δargument == Difference<Argument>{} &&
      values.first == values.second) {
    a2_ = {};
    return;
  }
  auto const one_over_Δargument = 1.0 / Δargument;
  auto const one_over_Δargument² = one_over_Δargument * one_over_Δargument;
  Difference<Value> const Δvalue = values.second - values.first;
  a2_ = Δvalue * one_over_Δargument² - derivative_first * one_over_Δargument;
}

template<typename Value, typename Argument>
Value Hermite2<Value, Argument>::Evaluate(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return ((a2_ * Δargument) + a1_) * Δargument + a0_;
}

template<typename Value, typename Argument>
typename Hermite2<Value, Argument>::Derivative1
Hermite2<Value, Argument>::EvaluateDerivative(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return 2.0 * a2_ * Δargument + a1_;
}

template<typename Value, typename Argument>
Argument Hermite2<Value, Argument>::FindExtremum() const {
  return arguments_.first - a1_ / (2 * a2_);
}

}  // namespace internal
}  // namespace _hermite2
}  // namespace numerics
}  // namespace principia
