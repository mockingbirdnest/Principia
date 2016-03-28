
#pragma once

#include "numerics/hermite3.hpp"

#include "numerics/root_finders.hpp"

namespace principia {
namespace numerics {

template<typename Argument, typename Value>
Hermite3<Argument, Value>::Hermite3(
    std::pair<Argument, Argument> const& arguments,
    std::pair<Value, Value> const& values,
    std::pair<Derivative, Derivative> const& derivatives)
    : arguments_(arguments) {
  a0_ = values.first;
  a1_ = derivatives.first;
  Difference<Argument> const Δargument = arguments_.second - arguments_.first;
  auto const one_over_Δargument = 1.0 / Δargument;
  auto const one_over_Δargument_squared =
      one_over_Δargument * one_over_Δargument;
  auto const one_over_Δargument_cubed =
      one_over_Δargument * one_over_Δargument_squared;
  Difference<Value> const Δvalue = values.second - values.first;
  a2_ = 3.0 * Δvalue * one_over_Δargument_squared -
            (2.0 * derivatives.first + derivatives.second) * one_over_Δargument;
  a3_ = -2.0 * Δvalue * one_over_Δargument_cubed +
            (derivatives.first + derivatives.second) *
                one_over_Δargument_squared;
}

template <typename Argument, typename Value>
Value Hermite3<Argument, Value>::Evaluate(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return (((a3_ * Δargument + a2_) * Δargument) + a1_) * Δargument + a0_;
}

template<typename Argument, typename Value>
typename Hermite3<Argument, Value>::Derivative
Hermite3<Argument, Value>::EvaluateDerivative(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return ((3.0 * a3_ * Δargument + 2.0 * a2_) * Δargument) + a1_;
}

template<typename Argument, typename Value>
std::set<Argument> Hermite3<Argument, Value>::FindExtrema() const {
  return SolveQuadraticEquation<Argument, Derivative>(
      arguments_.first, a1_, 2.0 * a2_, 3.0 * a3_);
}

}  // namespace numerics
}  // namespace principia
