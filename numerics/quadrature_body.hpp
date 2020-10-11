
#pragma once

#include "numerics/quadrature.hpp"

#include <memory>
#include <vector>

#include "mathematica/mathematica.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/gauss_legendre_weights.mathematica.h"
#include "numerics/legendre_roots.mathematica.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace quadrature {
namespace internal_quadrature {

inline mathematica::Logger logger(TEMP_DIR / "quadrature.wl",
                                  /*make_unique=*/false);

using base::FloorLog2;
using quantities::Angle;
using quantities::Cos;
using quantities::Difference;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

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

template<int points, typename Argument, typename Function, typename Function2>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtisImplementation(
    Function2 const& f2,
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    double const relative_tolerance,
    Primitive<std::invoke_result_t<Function, Argument>, Argument> const
        previous_estimate) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;
  Result const estimate = ClenshawCurtis<points>(f, lower_bound, upper_bound);
  // This is the naïve estimate mentioned in [Gen72b], p. 339.
  double const relative_error_estimate =
      Hilbert<Result>::Norm(previous_estimate - estimate) /
      Hilbert<Result>::Norm(estimate);
  // We look for an estimated relative error smaller than
  // |relative_tolerance * points|: since the integral is computed from |points|
  // evaluations of |f|, it will necessarily carry a relative error proportional
  // to |points|, so it makes no sense to look for convergence beyond that.
  logger.Append("estimate", std::tuple(points, estimate),
        mathematica::ExpressIn(Metre, Second, Radian));
  if (points > 1e6) {
    logger.Append(
        "function",
        std::tuple{f2, lower_bound, upper_bound, previous_estimate, estimate},
        mathematica::ExpressIn(Metre, Second, Radian));
  }
  if (relative_error_estimate > relative_tolerance * points) {
    if constexpr (points > 1 << 24) {
      LOG(FATAL) << "Too many refinements while integrating from "
                 << lower_bound << " to " << upper_bound;
    } else {
      return AutomaticClenshawCurtisImplementation<points + points - 1>(
          f2, f, lower_bound, upper_bound, relative_tolerance, estimate);
    }
  }
  return estimate;
}

template<int initial_points,
         typename Argument, typename Function, typename Function2>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtis(
    Function2 const& f2,
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    double const relative_tolerance) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;
  // TODO(egg): factor the evaluations of f.
  Result const estimate =
      ClenshawCurtis<initial_points>(f, lower_bound, upper_bound);
  return AutomaticClenshawCurtisImplementation<
      initial_points + initial_points - 1>(
      f2, f, lower_bound, upper_bound, relative_tolerance, estimate);
}

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound) {
  // We follow the notation from [Gen72b] and [Gen72c].
  using Value = std::invoke_result_t<Function, Argument>;

  constexpr int N = points - 1;
  constexpr int log2_N = FloorLog2(N);
  static_assert(N == 1 << log2_N);

  Difference<Argument> const half_width = (upper_bound - lower_bound) / 2;

  constexpr Angle N⁻¹π = π * Radian / N;

  // We use a discrete Fourier transform rather than a cosine transform, see
  // [Gen72c], equation (3).
  // TODO(egg): Consider a discrete cosine transform, ideally incrementally
  // computed to improve the performance of the automatic quadrature.

  // If we identify [lower_bound, upper_bound] with [-1, 1],
  // f_cos_N⁻¹π[s] is f(cos πs/N).
  std::vector<Value> f_cos_N⁻¹π;
  f_cos_N⁻¹π.resize(2 * N);
  double bigrel = 0;
  double bigs = -1;
  Value bignext;
  for (int s = 0; s <= N; ++s) {
    // The N + 1 evaluations.
    f_cos_N⁻¹π[s] = f(lower_bound + half_width * (1 + Cos(N⁻¹π * s)));
    auto const next = f(lower_bound + half_width * (1 + Cos(N⁻¹π * s) + 1.0e-7));
    if (f_cos_N⁻¹π[s] != Value{} && next != Value{}) {
      auto const rel = std::abs(1 - f_cos_N⁻¹π[s] / next);
      if (rel > bigrel) {
        bigrel = rel;
        bigs = s;
        bignext = next;
      }
    }
  }
  LOG(ERROR) << bigs << "/" << N << " " << bigrel << " " << f_cos_N⁻¹π[bigs]
             << " " << bignext << " "
             << lower_bound + half_width * (1 + Cos(N⁻¹π * bigs) + 1.0e-7);
  for (int s = N + 1; s <= 2 * N - 1; ++s) {
    f_cos_N⁻¹π[s] = f_cos_N⁻¹π[2 * N - s];  // [Gen72c] (5).
  }

  if (points > 1e6) {
    logger.Append("eval", f_cos_N⁻¹π,
                  mathematica::ExpressIn(Metre, Second, Radian));
  }

  auto const fft = std::make_unique<FastFourierTransform<Value, Angle, 2 * N>>(
      f_cos_N⁻¹π, N⁻¹π);
  auto const& a = *fft;

  // [Gen72b] equation (7), factoring out the division by N.
  Value Σʺ{};
  for (std::int64_t n = 0; n <= N; n += 2) {
    // The notation g is from [OLBC10], 3.5.17.
    int const gₙ = n == 0 || n == N ? 1 : 2;
    Σʺ += a[n].real_part() * gₙ / (1 - n * n);
  }
  Σʺ /= N;

  return Σʺ * half_width;
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
