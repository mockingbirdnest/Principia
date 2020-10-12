
#pragma once

#include "numerics/quadrature.hpp"

#include <memory>
#include <vector>

#include "base/bits.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/gauss_legendre_weights.mathematica.h"
#include "numerics/legendre_roots.mathematica.h"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace quadrature {
namespace internal_quadrature {

using base::BitReversedIncrement;
using base::FloorLog2;
using geometry::Hilbert;
using quantities::Angle;
using quantities::Cos;
using quantities::Difference;
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

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtisImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::optional<double> const max_relative_error,
    std::optional<int> const max_points,
    Primitive<std::invoke_result_t<Function, Argument>, Argument> const
        previous_estimate,
    std::vector<std::invoke_result_t<Function, Argument>>& cached_f_cos_N⁻¹π) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;

  Result const estimate =
      ClenshawCurtisImplementation<points>(
          f, lower_bound, upper_bound, cached_f_cos_N⁻¹π);

  // This is the naïve estimate mentioned in [Gen72b], p. 339.
  auto const absolute_error_estimate =
      Hilbert<Result>::Norm(previous_estimate - estimate);

  if ((!max_relative_error.has_value() ||
       absolute_error_estimate >
           max_relative_error.value() * Hilbert<Result>::Norm(estimate)) &&
      (!max_points.has_value() || points <= max_points.value())) {
    if constexpr (points > 1 << 24) {
      LOG(FATAL) << "Too many refinements while integrating from "
                 << lower_bound << " to " << upper_bound;
    } else {
      cached_f_cos_N⁻¹π.reserve(2 * points + 1);
      return AutomaticClenshawCurtisImplementation<2 * points>(
          f,
          lower_bound, upper_bound,
          max_relative_error, max_points,
          estimate,
          cached_f_cos_N⁻¹π);
    }
  }
  return estimate;
}

//TODO(phl):comment
template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
ClenshawCurtisImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::vector<std::invoke_result_t<Function, Argument>>& cached_f_cos_N⁻¹π) {
  // We follow the notation from [Gen72b] and [Gen72c].
  using Value = std::invoke_result_t<Function, Argument>;

  constexpr int N = points;
  constexpr int log2_N = FloorLog2(N);

  Difference<Argument> const half_width = (upper_bound - lower_bound) / 2;
  constexpr Angle N⁻¹π = π * Radian / N;

  FillClenshawCurtisCache<points>(
      f, lower_bound, upper_bound, cached_f_cos_N⁻¹π);

  //TODO(phl): See if the copy can be avoided.
  std::vector<Value> f_cos_N⁻¹π;
  f_cos_N⁻¹π.resize(2 * N);
  // An index in f_cos_N⁻¹π, corresponding to the increasing order of s for this
  // value of N.
  int reverse = 0;
  // An index in cached_f_cos_N⁻¹π, corresponding to the order in which the
  // values were put in the cache.  Note that entry 0 is special and holds the
  // value corresponding to s = N.
  for (int direct = 1; direct <= N; ++direct) {
    f_cos_N⁻¹π[reverse] = cached_f_cos_N⁻¹π[direct];
    if (reverse > 0) {
      f_cos_N⁻¹π[2 * N - reverse] = cached_f_cos_N⁻¹π[direct];  // [Gen72c] (5).
    }
    reverse = BitReversedIncrement(reverse, log2_N);
  }
  f_cos_N⁻¹π[N] = cached_f_cos_N⁻¹π[0];

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

template<int points, typename Argument, typename Function>
void FillClenshawCurtisCache(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::vector<std::invoke_result_t<Function, Argument>>& cached_f_cos_N⁻¹π) {
  //TODO(phl):Order
  // If we identify [lower_bound, upper_bound] with [-1, 1],
  // cached_f_cos_N⁻¹π contains f(cos πs/N).
  // We use a discrete Fourier transform rather than a cosine transform, see
  // [Gen72c], equation (3).

  //TODO(phl): points is a power of 2.
  constexpr int N = points;
  static_assert(N >= 1);
  constexpr int log2_N = FloorLog2(N);
  static_assert(N == 1 << log2_N);

  // The cache already has all the data we need.
  if (cached_f_cos_N⁻¹π.size() > N) {
    return;
  }

  Difference<Argument> const half_width = (upper_bound - lower_bound) / 2;
  constexpr Angle N⁻¹π = π * Radian / N;

  if constexpr (N == 1) {
    DCHECK(cached_f_cos_N⁻¹π.empty());
    //TODO(phl):comment this hack
    cached_f_cos_N⁻¹π.push_back(
        f(lower_bound + half_width * (1 + Cos(N⁻¹π * /*s=*/N))));
    cached_f_cos_N⁻¹π.push_back(
        f(lower_bound + half_width * (1 + Cos(N⁻¹π * /*s=*/0))));
  } else {
    FillClenshawCurtisCache<N / 2>(
        f, lower_bound, upper_bound, cached_f_cos_N⁻¹π);
    // N/2 evaluations for f(cos πs/N) with s odd: the values for s even have
    // already been computed by the previous recursive call.
    for (int s = 1; s < N; s += 2) {
      cached_f_cos_N⁻¹π.push_back(
          f(lower_bound + half_width * (1 + Cos(N⁻¹π * s))));
    }
  }
}

template<int initial_points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::optional<double> const max_relative_error,
    std::optional<int> const max_points) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;
  using Value = std::invoke_result_t<Function, Argument>;
  constexpr int N = initial_points;

  std::vector<Value> cached_f_cos_N⁻¹π;
  cached_f_cos_N⁻¹π.reserve(2 * N + 1);
  Result const estimate = ClenshawCurtisImplementation<N>(
      f, lower_bound, upper_bound, cached_f_cos_N⁻¹π);
  return AutomaticClenshawCurtisImplementation<2 * N>(
      f,
      lower_bound, upper_bound,
      max_relative_error, max_points,
      estimate,
      cached_f_cos_N⁻¹π);
}

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound) {
  using Value = std::invoke_result_t<Function, Argument>;
  constexpr int N = points;

  std::vector<Value> cached_f_cos_N⁻¹π;
  cached_f_cos_N⁻¹π.reserve(N + 1);
  return ClenshawCurtisImplementation<N>(
      f, lower_bound, upper_bound, cached_f_cos_N⁻¹π);
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
