#pragma once

#include "numerics/hermite3.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "numerics/root_finders.hpp"

namespace principia {
namespace numerics {
namespace _hermite3 {
namespace internal {

using namespace principia::numerics::_root_finders;

template<typename Value_, typename Argument_>
Hermite3<Value_, Argument_>::Hermite3(
    std::pair<Argument, Argument> arguments,
    std::pair<Value, Value> const& values,
    std::pair<Derivative1, Derivative1> const& derivatives)
    : arguments_(std::move(arguments)) {
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
  auto const one_over_Δargument² = one_over_Δargument * one_over_Δargument;
  auto const one_over_Δargument³ =
      one_over_Δargument * one_over_Δargument²;
  Difference<Value> const Δvalue = values.second - values.first;
  a2_ = 3.0 * Δvalue * one_over_Δargument² -
            (2.0 * derivatives.first + derivatives.second) * one_over_Δargument;
  a3_ = -2.0 * Δvalue * one_over_Δargument³ +
            (derivatives.first + derivatives.second) * one_over_Δargument²;
}

template<typename Value_, typename Argument_>
Value_ Hermite3<Value_, Argument_>::Evaluate(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return (((a3_ * Δargument + a2_) * Δargument) + a1_) * Δargument + a0_;
}

template<typename Value_, typename Argument_>
typename Hermite3<Value_, Argument_>::Derivative1
Hermite3<Value_, Argument_>::
EvaluateDerivative(Argument const& argument) const {
  Difference<Argument> const Δargument = argument - arguments_.first;
  return ((3.0 * a3_ * Δargument + 2.0 * a2_) * Δargument) + a1_;
}

template<typename Value_, typename Argument_>
BoundedArray<Argument_, 2> Hermite3<Value_, Argument_>::FindExtrema() const {
  return SolveQuadraticEquation<Argument, Derivative1>(
      arguments_.first, a1_, 2.0 * a2_, 3.0 * a3_);
}

template<typename Value_, typename Argument_>
template<typename Samples>
auto Hermite3<Value_, Argument_>::LInfinityError(
    Samples const& samples,
    std::function<Argument const&(typename Samples::value_type const&)> const&
        get_argument,
    std::function<Value const&(typename Samples::value_type const&)> const&
        get_value) const -> NormType {
  NormType result{};
  for (const auto& sample : samples) {
    result = std::max(result,
                      Hilbert<Difference<Value>>::Norm(
                          Evaluate(get_argument(sample)) - get_value(sample)));
  }
  return result;
}

template<typename Value_, typename Argument_>
template<typename Samples>
bool Hermite3<Value_, Argument_>::LInfinityErrorIsWithin(
    Samples const& samples,
    std::function<Argument const&(typename Samples::value_type const&)> const&
        get_argument,
    std::function<Value const&(typename Samples::value_type const&)> const&
        get_value,
    NormType const& tolerance) const {
  for (const auto& sample : samples) {
    if (Hilbert<Difference<Value>>::Norm(Evaluate(get_argument(sample)) -
                                         get_value(sample)) >= tolerance) {
      return false;
    }
  }
  return true;
}

}  // namespace internal
}  // namespace _hermite3
}  // namespace numerics
}  // namespace principia
