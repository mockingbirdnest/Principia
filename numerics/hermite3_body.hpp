#pragma once

#include "numerics/hermite3.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/root_finders.hpp"

namespace principia {
namespace numerics {
namespace _hermite3 {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_root_finders;

template<typename T>
struct Splitter;

template<typename Scalar>
struct Splitter<R3Element<Scalar>> {
  using Value = Scalar;

  std::vector<Scalar> Split(R3Element<Scalar> const& r3_element) {
    return {r3_element.x, r3_element.y, r3_element.z};
  }
};

template<typename Scalar, typename Frame, int rank>
  requires(rank == 2 || rank == 3)
struct Splitter<Multivector<Scalar, Frame, rank>> {
  using Value = Scalar;

  std::vector<Scalar> Split(
      Multivector<Scalar, Frame, rank> const& multivector) {
    return Splitter<R3Element<Scalar>>(multivector.coordinates());
  }
};

template<affine Value_, affine Argument_>
Hermite3<Value_, Argument_>::Hermite3(
    std::pair<Argument, Argument> const& arguments,
    std::pair<Value, Value> const& values,
    std::pair<Derivative1, Derivative1> const& derivatives)
    : p_(MakePolynomial(arguments, values, derivatives)),
      pʹ_(p_.Derivative()) {}

template<affine Value_, affine Argument_>
Value_ Hermite3<Value_, Argument_>::Evaluate(Argument const& argument) const {
  return p_(argument);
}

template<affine Value_, affine Argument_>
typename Hermite3<Value_, Argument_>::Derivative1
Hermite3<Value_, Argument_>::
EvaluateDerivative(Argument const& argument) const {
  return pʹ_(argument);
}

template<affine Value_, affine Argument_>
void Hermite3<Value_, Argument_>::EvaluateWithDerivative(
    Argument const& argument,
    Value& value,
    Derivative1& derivative) const {
  return p_.EvaluateWithDerivative(argument, value, derivative);
}

template<affine Value_, affine Argument_>
BoundedArray<Argument_, 2> Hermite3<Value_, Argument_>::FindExtrema() const {
  auto const& coefficients = pʹ_.coefficients();
  return SolveQuadraticEquation<Argument, Derivative1>(
      p_.origin(),
      std::get<0>(coefficients),
      std::get<1>(coefficients),
      std::get<2>(coefficients));
}

template<affine Value_, affine Argument_>
BoundedArray<Argument_, 2> Hermite3<Value_, Argument_>::FindExtrema(
    Argument const& lower,
    Argument const& upper) const {
  auto const extrema = FindExtrema();
  BoundedArray<Argument, 2> result;
  for (auto const& extremum : extrema) {
    if (lower <= extremum && extremum <= upper) {
      result.push_back(extremum);
    }
  }
  return result;
}

template<affine Value_, affine Argument_>
auto Hermite3<Value_, Argument_>::LInfinityL₁NormUpperBound(
    Argument const& lower,
    Argument const& upper) const -> NormType {
  CHECK_LE(lower, upper);

  // First split the coefficients of our polynomial by dimension.  This part is
  // *not* coordinate-free.
  auto const& coefficients = p_.coefficients();
  using S = Splitter<Value>;
  using P = PolynomialInMonomialBasis<S::Value, Argument, 3>;
  auto const split_a0 = S::Split(std::get<0>(coefficients));
  auto const split_a1 = S::Split(std::get<1>(coefficients));
  auto const split_a2 = S::Split(std::get<2>(coefficients));
  auto const split_a3 = S::Split(std::get<3>(coefficients));

  // Build a polynomial for each dimension.
  std::vector<P> split_polynomials;
  for (std::int64_t i = 0; i < split_a0.size(); ++i) {
    split_polynomials.emplace_back(typename P::Coefficients{
        split_a0[i], split_a1[i], split_a2[i], split_a3[i]});
  }

  NormType norm{};
  // Find the extrema of each polynomial.
  for (auto const& pᵢ : split_polynomials) {
    auto const extrema = pᵢ.FindExtrema(lower, upper);
    // For each extremum, evaluate all the polynomials and compute the sum of
    // their absolute value.  This is pessimistic (and is the reason why this
    // function only returns an upper bound).
    for (auto const extremum : extrema) {
      NormType sum{};
      for (auto const pⱼ : split_polynomials) {
        sum += Abs(pⱼ(extremum));
      }
      norm = std::max(norm, sum);
    }
  }

  return norm;
}

template<affine Value_, affine Argument_>
template<typename Samples>
auto Hermite3<Value_, Argument_>::LInfinityL₂Error(
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

template<affine Value_, affine Argument_>
template<typename Samples>
bool Hermite3<Value_, Argument_>::LInfinityL₂ErrorIsWithin(
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

template<affine Value_, affine Argument_>
FORCE_INLINE(inline)
PolynomialInMonomialBasis<Value_, Argument_, 3>
Hermite3<Value_, Argument_>::MakePolynomial(
    std::pair<Argument, Argument> const& arguments,
    std::pair<Value, Value> const& values,
    std::pair<Derivative1, Derivative1> const& derivatives) {
  using Derivative2 = Derivative<Derivative1, Argument>;
  using Derivative3 = Derivative<Derivative2, Argument>;

  Value const a0 = values.first;
  Derivative1 const a1 = derivatives.first;
  Derivative2 a2;
  Derivative3 a3;

  Difference<Argument> const Δargument = arguments.second - arguments.first;
  // If we were given the same point twice, there is a removable singularity.
  // Otherwise, if the arguments are the same but not the values or the
  // derivatives, we proceed to merrily NaN away as we should.
  if (Δargument != Difference<Argument>{} ||
      values.first != values.second ||
      derivatives.first != derivatives.second) {
    auto const one_over_Δargument = 1.0 / Δargument;
    auto const one_over_Δargument² = one_over_Δargument * one_over_Δargument;
    auto const one_over_Δargument³ = one_over_Δargument * one_over_Δargument²;
    Difference<Value> const Δvalue = values.second - values.first;
    a2 = 3.0 * Δvalue * one_over_Δargument² -
         (2.0 * derivatives.first + derivatives.second) * one_over_Δargument;
    a3 = -2.0 * Δvalue * one_over_Δargument³ +
         (derivatives.first + derivatives.second) * one_over_Δargument²;
  }
  return PolynomialInMonomialBasis<Value, Argument, 3>({a0, a1, a2, a3},
                                                       arguments.first);
}

}  // namespace internal
}  // namespace _hermite3
}  // namespace numerics
}  // namespace principia
