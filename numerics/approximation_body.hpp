#pragma once

#include "numerics/approximation.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "base/tags.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

// Compute the interpolation matrix and cache it in a static variable.
template<int N>
FixedMatrix<double, N + 1, N + 1> const& ЧебышёвInterpolationMatrix() {
  static FixedMatrix<double, N + 1, N + 1> const ℐ = []() {
    FixedMatrix<double, N + 1, N + 1> ℐ(uninitialized);
    for (std::int64_t j = 0; j <= N; ++j) {
      double const pⱼ = j == 0 || j == N ? 2 : 1;
      for (std::int64_t k = 0; k <= N; ++k) {
        double const pₖ = k == 0 || k == N ? 2 : 1;
        ℐ(j, k) = 2 * Cos(π * j * k * Radian / N) / (pⱼ * pₖ * N);
      }
    }
    return ℐ;
  }();
  return ℐ;
}

template<int N, int max_degree, typename Argument, typename Function>
ЧебышёвSeries<Value<Argument, Function>, Argument>
ЧебышёвPolynomialInterpolantImplementation(
    Function const& f,
    Argument const& a,
    Argument const& b,
    Difference<Value<Argument, Function>> const& max_error,
    FixedVector<Value<Argument, Function>, N / 2 + 1> const& previous_fₖ,
    FixedVector<Value<Argument, Function>, N / 2 + 1> const& previous_aⱼ,
    Difference<Value<Argument, Function>>* const error_estimate) {
  // This implementation follows [Boy13], section 4 and appendix A.
  auto const midpoint = Barycentre(std::pair{a, b}, std::pair{0.5, 0.5});

  auto чебышёв_lobato_point =
      [&a, &b, &midpoint](std::int64_t const k) -> Argument {
    return 0.5 * (b - a) * Cos(π * k * Radian / N) + midpoint;
  };

  FixedVector<Value<Argument, Function>, N + 1> fₖ(uninitialized);

  // Reuse the previous evaluations of |f|.
  for (std::int64_t k = 0; k <= N / 2; ++k) {
    fₖ[2 * k] = previous_fₖ[k];
  }

  // Evaluate |f| for the new points.
  for (std::int64_t k = 1; k < N; k += 2) {
    fₖ[k] = f(чебышёв_lobato_point(k));
  }

  // Compute the coefficients of the Чебышёв polynomial.
  FixedMatrix<double, N + 1, N + 1> const& ℐⱼₖ =
      ЧебышёвInterpolationMatrix<N>();
  auto const aⱼ = ℐⱼₖ * fₖ;

  // Compute an upper bound for the error, based on the previous and new
  // polynomials.
  Difference<Value<Argument, Function>> current_error_estimate{};
  for (std::int64_t j = 0; j <= N / 2; ++j) {
    current_error_estimate += Abs(previous_aⱼ[j] - aⱼ[j]);
  }
  for (std::int64_t j = N / 2 + 1; j <= N; ++j) {
    current_error_estimate += Abs(aⱼ[j]);
  }

  if constexpr (N <= max_degree) {
    if (current_error_estimate > max_error) {
      // Note that this recursive call overflows the stack when
      // max_degree >= 256.  We could allocate on the heap, but then we don't
      // care about very high degree polynomials.
      return ЧебышёвPolynomialInterpolantImplementation<2 * N, max_degree>(
          f, a, b, max_error, fₖ, aⱼ, error_estimate);
    }
  }

  // Unlike [Boy13], section 4, we return the polynomial of the lower degree
  // that is within the |max_error| bound (that is, the one of degree N / 2).
  // Even though we computed the polynomial of degree N, returning it would
  // impose an unnecessary cost on the client (e.g., more costly evaluation). If
  // a client wants a more precise approximation, they just need to give a
  // smaller |max_error|.
  if (error_estimate != nullptr) {
    *error_estimate = current_error_estimate;
  }
  std::vector<Value<Argument, Function>> coefficients;
  std::copy(previous_aⱼ.begin(), previous_aⱼ.end(),
            std::back_inserter(coefficients));
  return ЧебышёвSeries<Value<Argument, Function>, Argument>(
      coefficients, a, b);
}

template<int max_degree, typename Argument, typename Function>
ЧебышёвSeries<Value<Argument, Function>, Argument>
ЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    Difference<Value<Argument, Function>>* const error_estimate) {
  auto const& a = lower_bound;
  auto const& b = upper_bound;
  auto const f_a = f(a);
  auto const f_b = f(b);
  FixedVector<Value<Argument, Function>, 2> const fₖ({f_b, f_a});
  FixedVector<Value<Argument, Function>, 2> const aⱼ(
      {0.5 * (f_b + f_a), 0.5 * (f_b - f_a)});

  return ЧебышёвPolynomialInterpolantImplementation</*N=*/2, max_degree>(
      f, a, b, max_error, fₖ, aⱼ, error_estimate);
}

template<int max_degree, typename Argument, typename Function>
std::vector<ЧебышёвSeries<Value<Argument, Function>, Argument>>
AdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    Difference<Value<Argument, Function>>* error_estimate) {
  // Try to build an interpolation over the entire interval.
  Difference<Value<Argument, Function>> full_error_estimate;
  auto full_interpolant = ЧебышёвPolynomialInterpolant<max_degree>(
      f, lower_bound, upper_bound, max_error, &full_error_estimate);
  if (full_error_estimate <= max_error) {
    // If the interpolant over the entire interval is within the desired error
    // bound, return it.
    if (error_estimate != nullptr) {
      *error_estimate = full_error_estimate;
    }
    std::vector<ЧебышёвSeries<Value<Argument, Function>, Argument>>
        interpolants;
    interpolants.emplace_back(std::move(full_interpolant));
    return interpolants;
  } else {
    // If the interpolant over the entire interval is not within the desired
    // error bound, subdivide the interval.
    Difference<Value<Argument, Function>> upper_error_estimate;
    Difference<Value<Argument, Function>> lower_error_estimate;
    auto const midpoint =
        Barycentre(std::pair(lower_bound, upper_bound), std::pair(1, 1));
    auto lower_interpolants =
        AdaptiveЧебышёвPolynomialInterpolant<max_degree>(
            f, lower_bound, midpoint, max_error, &lower_error_estimate);
    auto upper_interpolants =
        AdaptiveЧебышёвPolynomialInterpolant<max_degree>(
            f, midpoint, upper_bound, max_error, &upper_error_estimate);
    std::vector<ЧебышёвSeries<Value<Argument, Function>, Argument>>
        all_interpolants;
    std::move(lower_interpolants.begin(),
              lower_interpolants.end(),
              std::back_inserter(all_interpolants));
    std::move(upper_interpolants.begin(),
              upper_interpolants.end(),
              std::back_inserter(all_interpolants));
    if (error_estimate != nullptr) {
      *error_estimate = std::max(lower_error_estimate, upper_error_estimate);
    }
    return all_interpolants;
  }
}

}  // namespace internal
}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
