#pragma once

#include "numerics/hermite3.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/concepts.hpp"

namespace principia {
namespace numerics {
namespace _hermite3 {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_root_finders;
using namespace principia::quantities::_concepts;

template<typename T>
struct Splitter;

template<quantity Q>
struct Splitter<Q> {
  static constexpr std::int64_t dimension = 1;
  using Value = Q;

  template<typename T>
  static std::array<T, dimension> Split(T const& t) {
    return {t};
  }
};

template<typename Scalar>
struct Splitter<R3Element<Scalar>> {
  static constexpr std::int64_t dimension = 3;
  using Value = Scalar;

  template<typename T>
  static std::array<T, dimension> Split(
      R3Element<T> const& r3_element) {
    return {r3_element.x, r3_element.y, r3_element.z};
  }
};

template<typename Scalar, typename Frame, int rank>
  requires(rank == 1 || rank == 2)
struct Splitter<Multivector<Scalar, Frame, rank>> {
  static constexpr std::int64_t dimension = 3;
  using Value = Scalar;

  template<typename T>
  static std::array<T, dimension> Split(
      Multivector<T, Frame, rank> const& multivector) {
    return Splitter<R3Element<T>>::Split(multivector.coordinates());
  }
};

template<affine Value_, affine Argument_>
Hermite3<Value_, Argument_>::Hermite3(
    std::pair<Argument, Argument> const& arguments,
    std::pair<Value, Value> const& values,
    std::pair<Derivative1, Derivative1> const& derivatives)
    : lower_(arguments.first),
      upper_(arguments.second),
      p_(MakePolynomial(arguments, values, derivatives)),
      pʹ_(p_.Derivative()) {}

template<affine Value_, affine Argument_>
Value_ Hermite3<Value_, Argument_>::operator()(Argument const& argument) const {
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
  auto const roots =
      SolveQuadraticEquation<Argument, Derivative1>(p_.origin(),
                                                    std::get<0>(coefficients),
                                                    std::get<1>(coefficients),
                                                    std::get<2>(coefficients));
  BoundedArray<Argument, 2> valid_roots;
  for (auto const& root : roots) {
    if (lower_ <= root && root <= upper_) {
      valid_roots.push_back(root);
    }
  }
  return valid_roots;
}

template<affine Value_, affine Argument_>
auto Hermite3<Value_, Argument_>::LInfinityL₁NormUpperBound() const
    -> NormType {
  // First split the coefficients of `p_` by dimension.  This part is *not*
  // coordinate-free.
  auto const& coefficients = p_.coefficients();
  auto const& origin = p_.origin();
  using S = Splitter<Value>;
  using H = Hermite3<typename S::Value, Argument>;
  using P = PolynomialInMonomialBasis<typename S::Value, Argument, degree>;

  // NOTE(phl): This could be done for any degree by shaving the tuple using
  // template metaprogramming, but our degree is 3, so unrolling is simpler.
  auto const split_a0 = S::Split(std::get<0>(coefficients));
  auto const split_a1 = S::Split(std::get<1>(coefficients));
  auto const split_a2 = S::Split(std::get<2>(coefficients));
  auto const split_a3 = S::Split(std::get<3>(coefficients));

  // Build a split Hermite polynomial for each dimension.
  std::vector<H> split_hermites;
  for (std::int64_t i = 0; i < S::dimension; ++i) {
    split_hermites.push_back(
        H(lower_,
          upper_,
          P({split_a0[i], split_a1[i], split_a2[i], split_a3[i]}, origin)));
  }

  // Find the extrema of each split polynomial.
  NormType sum{};
  for (H const& split_hermite : split_hermites) {
    auto const extrema = split_hermite.FindExtrema();
    // Compute the maximum of the absolute value of each split polynomial at its
    // extrema.  This is the L∞ norm of that polynomial.
    NormType max{};
    for (auto const extremum : extrema) {
      max = std::max(max, Abs(split_hermite(extremum)));
    }
    // The L₁ norm of `p_` is bounded by the sum of the L∞ norms of the split
    // polynomials.
    sum += max;
  }

  return sum;
}

template<affine Value_, affine Argument_>
auto Hermite3<Value_, Argument_>::LInfinityL₂Norm() const -> NormType {
  CHECK_EQ((*this)(lower_), Value{});
  CHECK_EQ(this->EvaluateDerivative(lower_), Derivative1{});

  // The value type of the polynomial `q(t) = p_(t) / (t - lower_)²`.
  using QValue = Quotient<Value, Square<Difference<Argument>>>;

  auto const& coefficients = p_.coefficients();
  PolynomialInMonomialBasis<QValue, Argument, 1> const q(
      {std::get<2>(coefficients), std::get<3>(coefficients)},
      lower_);
  // The monomial `(t - lower_)`.
  PolynomialInMonomialBasis<Difference<Argument>, Argument, 1> const monomial(
      {Difference<Argument>{}, 1.0}, lower_);

  // This is not quite `d/dt(‖p_(t)‖₂²)`, but it differs from it by a factor
  // `(t - lower)³`, which is irrelevant to find its zeroes.
  using DValue =
      Quotient<Square<NormType>, Exponentiation<Difference<Argument>, 4>>;
  PolynomialInMonomialBasis<DValue, Argument, 2> const norm₂²_derivative =
      PointwiseInnerProduct(q, (2.0 * q + monomial * q.Derivative()));
  auto const& norm₂²_derivative_coefficients = norm₂²_derivative.coefficients();
  auto const norm₂²_derivative_roots =
      SolveQuadraticEquation<Argument, DValue>(
          q.origin(),
          std::get<0>(norm₂²_derivative_coefficients),
          std::get<1>(norm₂²_derivative_coefficients),
          std::get<2>(norm₂²_derivative_coefficients));

  // The extrema of `‖p_(t)‖₂` are either at the roots of its derivative or at
  // the upper bound of the interval.
  Square<NormType> max = Hilbert<Difference<Value>>::Norm²((*this)(upper_));
  for (auto const& root : norm₂²_derivative_roots) {
    if (lower_ < root && root < upper_) {
      max = std::max(max, Hilbert<Difference<Value>>::Norm²((*this)(root)));
    }
  }

  return Sqrt(max);
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
                          (*this)(get_argument(sample)) - get_value(sample)));
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
    if (Hilbert<Difference<Value>>::Norm((*this)(get_argument(sample)) -
                                         get_value(sample)) >= tolerance) {
      return false;
    }
  }
  return true;
}

template<affine Value_, affine Argument_>
Hermite3<Value_, Argument_>::Hermite3(
    Argument const& lower,
    Argument const& upper,
    PolynomialInMonomialBasis<Value, Argument, degree> p)
    : lower_(lower),
      upper_(upper),
      p_(std::move(p)),
      pʹ_(p_.Derivative()) {
  CHECK_EQ(p_.origin(), lower_);
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
