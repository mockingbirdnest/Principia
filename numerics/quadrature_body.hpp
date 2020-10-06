
#pragma once

#include "numerics/quadrature.hpp"

#include "numerics/gauss_legendre_weights.mathematica.h"
#include "numerics/legendre_roots.mathematica.h"
#include "quantities/elementary_functions.hpp"
#include "numerics/double_precision.hpp"

namespace principia {
namespace numerics {
namespace quadrature {
namespace internal_quadrature {

using quantities::Abs;
using quantities::Cos;
using quantities::Difference;
using quantities::Pow;
using quantities::si::Radian;

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> Gauss(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    double const* const nodes,
    double const* const weights) {
  Difference<Argument> half_width = (upper_bound - lower_bound) / 2;
  std::invoke_result_t<Function, Argument> result{};
  for (int i = 0; i < points; ++i) {
    Argument const scaled_node = lower_bound + half_width * (nodes[i] + 1);
    // TODO(phl): Consider compensated summation.
    result += weights[i] * f(scaled_node);
  }
  return result * half_width;
}

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> GaussLegendre(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound) {
  static_assert(points < LegendreRoots.size,
                "No table for Gauss-Legendre with the chosen number of points");
  return Gauss<points>(f,
                       lower_bound,
                       upper_bound,
                       LegendreRoots[points],
                       GaussLegendreWeights[points]);
}

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    typename Hilbert<Primitive<std::invoke_result_t<Function, Argument>,
                               Argument>>::NormType const absolute_tolerance,
    double const relative_tolerance) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;
  // TODO(egg): factor the evaluations of f (and of cos).
  Result previous_estimate = ClenshawCurtis(f, lower_bound, upper_bound, 3);
  Result estimate = ClenshawCurtis(f, lower_bound, upper_bound, 5);
  std::int64_t last_points = 5;
  for (std::int64_t points = 9;
       Abs(previous_estimate / estimate - 1) > relative_tolerance &&
       Hilbert<Result>::Norm(previous_estimate - estimate) > absolute_tolerance;
       points += points - 1) {
    LOG(ERROR) << "..." << Abs(previous_estimate / estimate - 1) << ", "
               << previous_estimate - estimate
               << ". Clenshaw-Curtis with " << points << "...";
    previous_estimate = estimate;
    estimate = ClenshawCurtis(f, lower_bound, upper_bound, points);
    last_points = points;
  }
  LOG(ERROR) << "Clenshaw-Curtis used " << last_points
             << " points for a relative error of "
             << Abs(previous_estimate / estimate - 1) << ", absolute error of "
             << previous_estimate - estimate;
  CHECK_EQ(estimate, estimate);
  return estimate;
}

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::int64_t const points) {
  std::int64_t const n = points - 1;
  Difference<Argument> const half_width = (upper_bound - lower_bound) / 2;
  // This vector contains nodes, which are the extrema of the Чебышёв polynomial
  // Tn, together with ±1: cos_n⁻¹π[k] = cos kπ/n.
  // TODO(egg): Consider caching this vector, we will keep coming here with n
  // equal to a power of two (we should probably share it with the FFT as well).
  std::vector<double> cos_n⁻¹π;
  cos_n⁻¹π.resize(n + 1);
  for (std::int64_t k = 0; k <= n; ++k) {
    cos_n⁻¹π[k] = Cos(k * π * Radian / n);
  }

  // See [OLBC10], 3.5(iv), and in particular 3.5.16 and 3.5.17 therein.
  DoublePrecision<std::invoke_result_t<Function, Argument>> Σₖ_wₖ_f_xₖ{};
  for (std::int64_t k = 0; k <= n; ++k) {
    Argument const xₖ = lower_bound + half_width * (cos_n⁻¹π[k] + 1);

    double const gₖ = k == 0 || k == n ? 1 : 2;
    DoublePrecision<double> Σⱼ;
    for (std::int64_t j = 1; j <= n / 2; ++j) {
      double const bⱼ = 2 * j == n ? 1 : 2;
      std::int64_t const jk_mod_n = (j * k) % n;
      double const cos_2jkπn⁻¹ = 2 * jk_mod_n <= n
                                     ? cos_n⁻¹π[2 * jk_mod_n]
                                     : cos_n⁻¹π[2 * (n - jk_mod_n)];
      Σⱼ += DoublePrecision<double>((bⱼ / (4 * (j * j) - 1) * cos_2jkπn⁻¹));
    }
    double const wₖ = gₖ / n * (1 - Σⱼ.value);

    Σₖ_wₖ_f_xₖ += DoublePrecision(wₖ * f(xₖ));
  }

  return Σₖ_wₖ_f_xₖ.value * half_width;
}

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> Midpoint(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    int const intervals) {
  Difference<Argument> const h = (upper_bound - lower_bound) / intervals;
  Primitive<std::invoke_result_t<Function, Argument>, Argument> result{};
  for (int i = 0; i < intervals; ++i) {
    result += f(lower_bound + (i + 0.5) * h) * h;
  }
  return result;
}

}  // namespace internal_quadrature
}  // namespace quadrature
}  // namespace numerics
}  // namespace principia
