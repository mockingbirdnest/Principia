
#pragma once

#include "numerics/hermite3.hpp"

#include <utility>
#include <vector>

#include "numerics/root_finders.hpp"

namespace principia {
namespace numerics {
namespace internal_hermite3 {

using quantities::Difference;

template<typename Argument, typename Value>
Hermite3<Argument, Value>::Hermite3(
    std::pair<Argument, Argument> const& arguments,
    std::pair<Value, Value> const& values,
    std::pair<Derivative1, Derivative1> const& derivatives)
    : arguments_(arguments) {
  a0_ = values.first;
  a1_ = derivatives.first;
  Difference<Argument> const Δargument = arguments_.second - arguments_.first;
  // If we were given the same point twice, there is a removable singularity.
  // Otherwise, if the arguments are the same but not the values or the
  // derivatives, we proceed to merrily NaN away as we should.
  if (Δargument == Difference<Argument>{} &&
      values.first == values.second &&
      derivatives.first == derivatives.second) {
    a2_ = {};
    a3_ = {};
    return;
  }
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

template<typename Argument, typename Value>
Value Hermite3<Argument, Value>::Evaluate(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return (((a3_ * Δargument + a2_) * Δargument) + a1_) * Δargument + a0_;
}

template<typename Argument, typename Value>
typename Hermite3<Argument, Value>::Derivative1
Hermite3<Argument, Value>::EvaluateDerivative(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return ((3.0 * a3_ * Δargument + 2.0 * a2_) * Δargument) + a1_;
}

template<typename Argument, typename Value>
BoundedArray<Argument, 2> Hermite3<Argument, Value>::FindExtrema() const {
  return SolveQuadraticEquation<Argument, Derivative1>(
      arguments_.first, a1_, 2.0 * a2_, 3.0 * a3_);
}

}  // namespace internal_hermite3
}  // namespace numerics
}  // namespace principia
