#pragma once

#include "numerics/quadrature.hpp"

#include <algorithm>
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

using namespace principia::base::_bits;
using namespace principia::geometry::_hilbert;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

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
void FillClenshawCurtisCache(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::vector<std::invoke_result_t<Function, Argument>>&
        f_cos_N⁻¹π_bit_reversed);

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
    std::vector<std::invoke_result_t<Function, Argument>>&
        f_cos_N⁻¹π_bit_reversed);

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
ClenshawCurtisImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::vector<std::invoke_result_t<Function, Argument>>&
        f_cos_N⁻¹π_bit_reversed);

// Our automatic Cleshaw-Curtis implementation doubles the number of points
// repeatedly until it reaches a suitable exit criterion.  Naïvely evaluating
// the function N times for each iteration would be wasteful.  Assume that we
// have already evaluated the function at 0/2 and 1/2 (that is, s = 0 and 1 for
// N = 2), we would evaluate it again at 0/4 and 2/4 for N = 4, and then proceed
// with evaluations at 1/4 and 3/4.  Overall, we would do roughly twice the
// work.
// Therefore, we cache the past evaluations.  However, we do not want to use a
// map and we do not want to do insertion in a vector, say, to squeeze 1/4
// between 0/2 and 1/2.  So we append the result of evaluations to the cache
// vector as we perform them.  To fill the cache for N points, we recursively
// fill it for N/2 points (which yields all the values at s/N for s even) and
// then we append the values at s/N for s odd.  The order of the second part
// matters.  It is done by taking the indices of the first part (covering
// [0/N, (N/2-1)/N[) in bit-reversed order and for each index s appending the
// value at (2 s + 1)/N.  This ensures that the entire array follows
// bit-reversed ordering.
// There is an extra wart: we need the value at N/N, but it does not fit nicely
// in the bit-reversed ordering.  So we store it at position 0 in the cache, and
// the regular cache starts at index 1.
// Clients are expected to reserve (at least) points entries in the cache vector
// for efficient heap allocation.
template<int points, typename Argument, typename Function>
void FillClenshawCurtisCache(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::vector<std::invoke_result_t<Function, Argument>>&
        f_cos_N⁻¹π_bit_reversed) {
  // If we identify [lower_bound, upper_bound] with [-1, 1],
  // f_cos_N⁻¹π_bit_reversed contains f(cos πs/N) in bit-reversed order of s.
  // We use a discrete Fourier transform rather than a cosine transform, see
  // [Gen72c], equation (3).

  // The cache already has all the data we need.
  if (f_cos_N⁻¹π_bit_reversed.size() >= points) {
    return;
  }

  constexpr int N = points - 1;
  static_assert(N >= 1);
  constexpr int log2_N = FloorLog2(N);
  static_assert(N == 1 << log2_N);

  Difference<Argument> const half_width = (upper_bound - lower_bound) / 2;
  constexpr Angle N⁻¹π = π * Radian / N;

  if constexpr (N == 1) {
    DCHECK(f_cos_N⁻¹π_bit_reversed.empty());
    // See above for the magic entry at index 0.  The entry at index 1 is a bona
    // fide entry which will be used by the recursive calls.
    f_cos_N⁻¹π_bit_reversed.push_back(f(lower_bound));  // s = N.
    f_cos_N⁻¹π_bit_reversed.push_back(f(upper_bound));  // s = 0.
  } else {
    // Fill the first half of the cache, corresponding to s even for this value
    // of N.
    FillClenshawCurtisCache<N / 2 + 1>(f,
                                       lower_bound, upper_bound,
                                       f_cos_N⁻¹π_bit_reversed);
    // N/2 evaluations for f(cos πs/N) with s odd.  Note the need to preserve
    // bit-reversed ordering.
    int reverse = 0;
    for (int evaluations = 0;
         evaluations < N / 2;
         ++evaluations, reverse = BitReversedIncrement(reverse, log2_N - 1)) {
      int const s = 2 * reverse + 1;
      f_cos_N⁻¹π_bit_reversed.push_back(
          f(lower_bound + half_width * (1 + Cos(N⁻¹π * s))));
    }
  }
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
    std::vector<std::invoke_result_t<Function, Argument>>&
        f_cos_N⁻¹π_bit_reversed) {
  using Result = Primitive<std::invoke_result_t<Function, Argument>, Argument>;

  Result const estimate =
      ClenshawCurtisImplementation<points>(
          f, lower_bound, upper_bound, f_cos_N⁻¹π_bit_reversed);

  // This is the naïve estimate mentioned in [Gen72b], p. 339.
  auto const absolute_error_estimate =
      Hilbert<Result>::Norm(previous_estimate - estimate);

  if ((!max_relative_error.has_value() ||
       absolute_error_estimate >
           max_relative_error.value() * Hilbert<Result>::Norm(estimate)) &&
      (!max_points.has_value() || points < max_points.value())) {
    if constexpr (points > 1 << 24) {
      LOG(FATAL) << "Too many refinements while integrating from "
                 << lower_bound << " to " << upper_bound
                 << ", relative error is "
                 << absolute_error_estimate / Hilbert<Result>::Norm(estimate);
    } else {
      f_cos_N⁻¹π_bit_reversed.reserve(2 * points - 1);
      return AutomaticClenshawCurtisImplementation<2 * points - 1>(
          f,
          lower_bound, upper_bound,
          max_relative_error, max_points,
          estimate,
          f_cos_N⁻¹π_bit_reversed);
    }
  }
  return estimate;
}

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
ClenshawCurtisImplementation(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::vector<std::invoke_result_t<Function, Argument>>&
        f_cos_N⁻¹π_bit_reversed) {
  // We follow the notation from [Gen72b] and [Gen72c].
  using Value = std::invoke_result_t<Function, Argument>;

  constexpr int N = points - 1;
  constexpr int log2_N = FloorLog2(N);

  Difference<Argument> const half_width = (upper_bound - lower_bound) / 2;
  constexpr Angle N⁻¹π = π * Radian / N;

  FillClenshawCurtisCache<points>(
      f, lower_bound, upper_bound, f_cos_N⁻¹π_bit_reversed);

  // TODO(phl): If might be possible to avoid copies since
  // f_cos_N⁻¹π_bit_reversed is tantalizing close to the order needed for the
  // FFT.
  std::vector<Value> f_cos_N⁻¹π;
  f_cos_N⁻¹π.resize(2 * N);
  // An index in f_cos_N⁻¹π in the order described by [Gen72b] and [Gen72c].
  int s = 0;
  // An index in f_cos_N⁻¹π_bit_reversed, corresponding to the order in which
  // the values were put in the cache.  Note that entry 0 is special and holds
  // the value corresponding to s = N.
  for (int bit_reversed_s = 1;
       bit_reversed_s <= N;
       ++bit_reversed_s, s = BitReversedIncrement(s, log2_N)) {
    f_cos_N⁻¹π[s] = f_cos_N⁻¹π_bit_reversed[bit_reversed_s];
    if (s > 0) {
      f_cos_N⁻¹π[2 * N - s] = f_cos_N⁻¹π[s];  // [Gen72c] (5).
    }
  }
  f_cos_N⁻¹π[N] = f_cos_N⁻¹π_bit_reversed[0];

  // TODO(phl): We could save some time by implementing a proper cosine
  // transform.
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
Primitive<std::invoke_result_t<Function, Argument>, Argument> GaussLegendre(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound) {
  static_assert(points < LegendreRoots.rows(),
                "No table for Gauss-Legendre with the chosen number of points");
  return Gauss<points>(f,
                       lower_bound,
                       upper_bound,
                       LegendreRoots.row<points>(),
                       GaussLegendreWeights.row<points>());
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
  std::vector<Value> f_cos_N⁻¹π_bit_reversed;
  f_cos_N⁻¹π_bit_reversed.reserve(2 * initial_points - 1);
  Result const estimate = ClenshawCurtisImplementation<initial_points>(
      f, lower_bound, upper_bound, f_cos_N⁻¹π_bit_reversed);
  return AutomaticClenshawCurtisImplementation<2 * initial_points - 1>(
      f,
      lower_bound, upper_bound,
      max_relative_error, max_points,
      estimate,
      f_cos_N⁻¹π_bit_reversed);
}

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound) {
  using Value = std::invoke_result_t<Function, Argument>;
  std::vector<Value> f_cos_N⁻¹π_bit_reversed;
  f_cos_N⁻¹π_bit_reversed.reserve(points);
  return ClenshawCurtisImplementation<points>(
      f, lower_bound, upper_bound, f_cos_N⁻¹π_bit_reversed);
}

inline std::optional<int> MaxPointsHeuristicsForAutomaticClenshawCurtis(
    AngularFrequency const& max_ω,
    Time const& Δt,
    int min_points_overall,
    int points_per_period) {
  return max_ω == AngularFrequency()
             ? std::optional<int>{}
             : std::max(min_points_overall,
                        static_cast<int>(points_per_period * Δt * max_ω /
                                         (2 * π * Radian)));
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
