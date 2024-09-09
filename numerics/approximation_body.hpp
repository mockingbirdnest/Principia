#pragma once

#include "numerics/approximation.hpp"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "base/tags.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/чебышёв_lobatto.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_чебышёв_lobatto;
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
not_null<std::unique_ptr<
    PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument>>>
ЧебышёвPolynomialInterpolantImplementation(
    Function const& f,
    Argument const& a,
    Argument const& b,
    Difference<Value<Argument, Function>> const& max_error,
    FixedVector<Value<Argument, Function>, N / 2 + 1> const& previous_fₖ,
    FixedVector<Value<Argument, Function>, N / 2 + 1> const& previous_aⱼ,
    Difference<Value<Argument, Function>>* const error_estimate) {
  // This implementation follows [Boy13], section 4 and appendix A.
  auto const midpoint = Barycentre({a, b});

  auto чебышёв_lobato_point =
      [&a, &b, &midpoint](std::int64_t const k) -> Argument {
    return 0.5 * (b - a) * ЧебышёвLobattoPoint<N>(k) + midpoint;
  };

  FixedVector<Value<Argument, Function>, N + 1> fₖ(uninitialized);

  // Reuse the previous evaluations of `f`.
  for (std::int64_t k = 0; k <= N / 2; ++k) {
    fₖ[2 * k] = previous_fₖ[k];
  }

  // Evaluate `f` for the new points.
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
  // that is within the `max_error` bound (that is, the one of degree N / 2).
  // Even though we computed the polynomial of degree N, returning it would
  // impose an unnecessary cost on the client (e.g., more costly evaluation). If
  // a client wants a more precise approximation, they just need to give a
  // smaller `max_error`.
  if (error_estimate != nullptr) {
    *error_estimate = current_error_estimate;
  }
  using Interpolant =
      PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument, N / 2>;
  return make_not_null_unique<Interpolant>(
      typename Interpolant::Coefficients(previous_aⱼ), a, b);
}

// Returns true if production of interpolants should stop.
template<int max_degree, typename Argument, typename Function>
bool StreamingAdaptiveЧебышёвPolynomialInterpolantImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    TerminationPredicate<Value<Argument, Function>, Argument> const& stop,
    Difference<Value<Argument, Function>>* const error_estimate) {
  // Try to build an interpolation over the entire interval.
  Difference<Value<Argument, Function>> full_error_estimate;
  auto full_interpolant = ЧебышёвPolynomialInterpolant<max_degree>(
      f, lower_bound, upper_bound, max_error, &full_error_estimate);
  if (full_error_estimate <= max_error ||
      !subdivide(*full_interpolant, full_error_estimate)) {
    // If the interpolant over the entire interval is within the desired error
    // bound, return it.  Same thing if `subdivide` tells us that we should not
    // subdivide the interval.
    VLOG(1) << "Degree " << full_interpolant->degree() << " interpolant over ["
            << lower_bound << " (" << f(lower_bound) << "), " << upper_bound
            << " (" << f(upper_bound) << ")] has error " << full_error_estimate;
    if (error_estimate != nullptr) {
      *error_estimate = full_error_estimate;
    }
    return stop(std::move(full_interpolant));
  } else {
    // If the interpolant over the entire interval is not within the desired
    // error bound, subdivide the interval.
    VLOG(1) << "Splitting [" << lower_bound << " (" << f(lower_bound) << "), "
            << upper_bound << " (" << f(upper_bound) << ")] with error "
            << full_error_estimate;
    Difference<Value<Argument, Function>> upper_error_estimate;
    Difference<Value<Argument, Function>> lower_error_estimate;
    auto const midpoint = Barycentre({lower_bound, upper_bound});
    bool const lower_interpolants_stop =
        StreamingAdaptiveЧебышёвPolynomialInterpolantImplementation<max_degree>(
            f,
            lower_bound,
            midpoint,
            max_error,
            subdivide,
            stop,
            &lower_error_estimate);
    if (lower_interpolants_stop) {
      if (error_estimate != nullptr) {
        *error_estimate = lower_error_estimate;
      }
      return true;
    }
    bool const upper_interpolants_stop =
        StreamingAdaptiveЧебышёвPolynomialInterpolantImplementation<max_degree>(
            f,
            midpoint,
            upper_bound,
            max_error,
            subdivide,
            stop,
            &upper_error_estimate);
    if (error_estimate != nullptr) {
      *error_estimate = std::max(lower_error_estimate, upper_error_estimate);
    }
    return upper_interpolants_stop;
  }
}

template<int max_degree, typename Argument, typename Function>
not_null<std::unique_ptr<
    PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument>>>
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
std::vector<not_null<std::unique_ptr<
    PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument>>>>
AdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    Difference<Value<Argument, Function>>* const error_estimate) {
  using Interpolant =
      PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument>;
  std::vector<not_null<std::unique_ptr<Interpolant>>> interpolants;

  TerminationPredicate<Value<Argument, Function>,
                       Argument> const emplace_back_and_continue =
      [&interpolants](not_null<std::unique_ptr<Interpolant>> interpolant) {
        interpolants.emplace_back(std::move(interpolant));
        return false;
      };

  StreamingAdaptiveЧебышёвPolynomialInterpolantImplementation<max_degree>(
      f,
      lower_bound,
      upper_bound,
      max_error,
      subdivide,
      emplace_back_and_continue,
      error_estimate);

  return interpolants;
}

template<int max_degree, typename Argument, typename Function>
void StreamingAdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    TerminationPredicate<Value<Argument, Function>, Argument> const& stop,
    Difference<Value<Argument, Function>>* const error_estimate) {
  StreamingAdaptiveЧебышёвPolynomialInterpolantImplementation<max_degree>(
      f,
      lower_bound,
      upper_bound,
      max_error,
      subdivide,
      stop,
      error_estimate);
}

}  // namespace internal
}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
