#pragma once

#include "numerics/approximation.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

constexpr std::int64_t max_чебышёв_degree = 500;

template<int N, typename Argument, typename Function>
ЧебышёвSeries<Value<Argument, Function>, Argument>
ЧебышёвPolynomialInterpolantImplementation(
    Function const& f,
    Argument const& a,
    Argument const& b,
    Difference<Value<Argument, Function>> const& max_error,
    FixedVector<Value<Argument, Function>, N / 2 + 1> const& previous_fₖ,
    FixedVector<Value<Argument, Function>, N / 2 + 1> const& previous_aⱼ) {
  auto const midpoint = Barycentre(std::pair{a, b}, std::pair{0.5, 0.5});

  auto чебышёв_lobato_point =
      [&a, &b, &midpoint](std::int64_t const k) -> Argument {
    return 0.5 * (b - a) * Cos(π * k * Radian / N) + midpoint;
  };

  FixedVector<Value<Argument, Function>, N + 1> fₖ;

  // Reuse the previous evaluations of |f|.
  for (std::int64_t k = 0; k <= N / 2; ++k) {
    fₖ[2 * k] = previous_fₖ[k];
  }

  // Evaluate |f| for the new points.
  for (std::int64_t k = 1; k < N; k += 2) {
    fₖ[k] = f(чебышёв_lobato_point(k));
  }

  FixedMatrix<double, N + 1, N + 1> ℐⱼₖ;
  for (std::int64_t j = 0; j <= N; ++j) {
    double const pⱼ = j == 0 || j == N ? 2 : 1;
    for (std::int64_t k = 0; k <= N; ++k) {
      double const pₖ = k == 0 || k == N ? 2 : 1;
      ℐⱼₖ(j, k) = 2 * Cos(π * j * k * Radian / N) / (pⱼ * pₖ * N);
    }
  }

  // Compute the coefficients of the Чебышёв polynomial.
  auto const aⱼ = ℐⱼₖ * fₖ;

  // Compute an upper bound for the error, based on the previous and new
  // polynomials.
  Difference<Value<Argument, Function>> error_estimate{};
  for (std::int64_t j = 0; j <= N / 2; ++j) {
    error_estimate += Abs(previous_aⱼ[j] - aⱼ[j]);
  }
  for (std::int64_t j = N / 2 + 1; j <= N; ++j) {
    error_estimate -= aⱼ[j];
  }

  LOG(ERROR)<<N<<" "<<error_estimate;
  if constexpr (2 * N <= max_чебышёв_degree) {
    if (error_estimate > max_error) {
      return ЧебышёвPolynomialInterpolantImplementation<2 * N>(
          f, a, b, max_error, fₖ, aⱼ);
    }
  }
  std::vector<Value<Argument, Function>> coefficients;
  std::copy(aⱼ.begin(), aⱼ.end(), std::back_inserter(coefficients));
  return ЧебышёвSeries<Value<Argument, Function>, Argument>(
      coefficients, a, b);
}

template<typename Argument, typename Function>
ЧебышёвSeries<Value<Argument, Function>, Argument>
ЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error) {
  auto const& a = lower_bound;
  auto const& b = upper_bound;
  auto const f_a = f(a);
  auto const f_b = f(b);
  FixedVector<Value<Argument, Function>, 2> const fₖ({f_b, f_a});
  FixedVector<Value<Argument, Function>, 2> const aⱼ(
      {0.5 * (f_b + f_a), 0.5 * (f_b - f_a)});
  return ЧебышёвPolynomialInterpolantImplementation</*N=*/2,
                                                    Argument,
                                                    Function>(
      f, a, b, max_error, fₖ, aⱼ);
}

}  // namespace internal
}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
