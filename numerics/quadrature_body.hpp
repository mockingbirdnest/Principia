
#pragma once

#include "numerics/quadrature.hpp"

#include <memory>
#include <vector>

#include "geometry/hilbert.hpp"
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
using geometry::Hilbert;
using quantities::Angle;
using quantities::Cos;
using quantities::Difference;
using quantities::Infinity;
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

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtisImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    double const relative_tolerance,
    Primitive<std::invoke_result_t<Function, Argument>, Argument> const
        previous_estimate) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;
  Result estimate;
  typename Hilbert<std::invoke_result_t<Function, Argument>>::NormType max;
  ClenshawCurtisImplementation<points>(f,
                                       lower_bound, upper_bound,
                                       estimate,
                                       max);

  // This is the naïve estimate mentioned in [Gen72b], p. 339.
  auto const absolute_error_estimate =
      Hilbert<Result>::Norm(previous_estimate - estimate);

  logger.Append("estimates",
                std::tuple(points, previous_estimate, estimate),
                mathematica::ExpressIn(Metre, Second, Radian));
  LOG(ERROR)<<points<<" "<<estimate<<" "<<previous_estimate<<" "
    <<absolute_error_estimate<<" "<<max<<" "<<
    relative_tolerance * max * (upper_bound - lower_bound);

  //TODO(phl):Comment
  // We look for an estimated relative error smaller than
  // |relative_tolerance * points|: since the integral is computed from |points|
  // evaluations of |f|, it will necessarily carry a relative error proportional
  // to |points|, so it makes no sense to look for convergence beyond that.
  if (absolute_error_estimate >
      relative_tolerance * max * (upper_bound - lower_bound)) {
    if constexpr (points > 1 << 24) {
      LOG(FATAL) << "Too many refinements while integrating from "
                 << lower_bound << " to " << upper_bound;
    } else {
      return AutomaticClenshawCurtisImplementation<points + points - 1>(
          f, lower_bound, upper_bound, relative_tolerance, estimate);
    }
  }
  return estimate;
}

template<int points, typename Argument, typename Function>
void ClenshawCurtisImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Primitive<std::invoke_result_t<Function, Argument>, Argument>& integral,
    typename Hilbert<std::invoke_result_t<Function, Argument>>::NormType& max) {
  // We follow the notation from [Gen72b] and [Gen72c].
  using Value = std::invoke_result_t<Function, Argument>;

  constexpr int N = points - 1;
  constexpr int log2_N = FloorLog2(N);
  static_assert(N == 1 << log2_N);

  max = -Infinity<typename Hilbert<Value>::NormType>;
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
  for (int s = 0; s <= N; ++s) {
    // The N + 1 evaluations.
    f_cos_N⁻¹π[s] = f(lower_bound + half_width * (1 + Cos(N⁻¹π * s)));
    max = std::max(max, Hilbert<Value>::Norm(f_cos_N⁻¹π[s]));
  }
  for (int s = N + 1; s <= 2 * N - 1; ++s) {
    f_cos_N⁻¹π[s] = f_cos_N⁻¹π[2 * N - s];  // [Gen72c] (5).
  }

  logger.Append("evaluations", std::tuple(points, f_cos_N⁻¹π),
                mathematica::ExpressIn(Metre, Second, Radian));

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

  integral = Σʺ * half_width;
}

template<int initial_points,
         typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtis(
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
      f, lower_bound, upper_bound, relative_tolerance, estimate);
}

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound) {
  Primitive<std::invoke_result_t<Function, Argument>, Argument> integral;
  typename Hilbert<std::invoke_result_t<Function, Argument>>::NormType max;
  ClenshawCurtisImplementation<points>(f,
                                       lower_bound, upper_bound,
                                       integral,
                                       max);
  return integral;
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
