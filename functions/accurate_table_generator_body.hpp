#pragma once

#include "functions/accurate_table_generator.hpp"

#include <algorithm>
#include <chrono>
#include <concepts>
#include <future>
#include <limits>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/strings/str_cat.h"
#include "base/bits.hpp"
#include "base/for_all_of.hpp"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.
#include "geometry/interval.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/lattices.hpp"
#include "numerics/matrix_views.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

using namespace principia::base::_bits;
using namespace principia::base::_for_all_of;
using namespace principia::geometry::_interval;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_lattices;
using namespace principia::numerics::_matrix_views;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;

constexpr std::int64_t T_max = 16;
static_assert(T_max >= 1);

template<std::int64_t zeroes>
bool HasDesiredZeroes(cpp_bin_float_50 const& y) {
  std::int64_t y_exponent;
  auto const y_mantissa = frexp(y, &y_exponent);
  auto const y_mantissa_scaled =
      ldexp(y_mantissa, std::numeric_limits<double>::digits);
  auto const y_post_mantissa = y_mantissa_scaled - floor(y_mantissa_scaled);
  auto const y_candidate_zeroes = ldexp(y_post_mantissa, zeroes);
  return trunc(y_candidate_zeroes) == 0;
}

template<std::int64_t zeroes, typename Container>
  requires std::same_as<typename Container::value_type, AccurateFunction>
bool AllFunctionValuesHaveDesiredZeroes(
    Container const& functions,
    cpp_rational const& argument) {
  return std::all_of(functions.begin(),
                     functions.end(),
                     [&argument](AccurateFunction const& f) {
                       return HasDesiredZeroes<zeroes>(f(argument));
                     });
}

struct StehléZimmermannSpecification {
  std::array<AccurateFunction, 2> functions;
  std::array<AccuratePolynomial<cpp_rational, 2>, 2> polynomials;
  std::array<AccurateFunction, 2> remainders;
  cpp_rational argument;
};

// Scales the argument, functions, polynomials, and remainders to lie within
// [1/2, 1[.
StehléZimmermannSpecification ScaleToBinade0(
    StehléZimmermannSpecification const& input,
    double& argument_scale) {
  auto const& functions = input.functions;
  auto const& polynomials = input.polynomials;
  auto const& remainders = input.remainders;
  auto const& starting_argument = input.argument;

  // TODO(phl): Handle an argument that is exactly a power of 2.
  std::int64_t argument_exponent;
  auto const argument_mantissa = frexp(
      static_cast<cpp_bin_float_50>(starting_argument), &argument_exponent);
  CHECK_NE(0, argument_mantissa);
  argument_scale = exp2(-argument_exponent);
  cpp_rational const scaled_argument = starting_argument * argument_scale;

  std::array<double, 2> function_scales;
  std::array<AccurateFunction, 2> scaled_functions;
  std::array<AccurateFunction, 2> scaled_remainders;
  for (std::int64_t i = 0; i < scaled_functions.size(); ++i) {
    std::int64_t function_exponent;
    auto const function_mantissa =
        frexp(functions[i](starting_argument), &function_exponent);
    CHECK_NE(0, function_mantissa);
    function_scales[i] = exp2(-function_exponent);
    scaled_functions[i] = [argument_scale,
                           function_scale = function_scales[i],
                           function = functions[i],
                           i](cpp_rational const& argument) {
      return function_scale * function(argument / argument_scale);
    };
    scaled_remainders[i] = [argument_scale,
                            function_scale = function_scales[i],
                            remainder = remainders[i],
                            i](cpp_rational const& argument) {
      return function_scale * remainder(argument / argument_scale);
    };
  }

  auto build_scaled_polynomial =
      [argument_scale, &starting_argument](
          double const function_scale,
          AccuratePolynomial<cpp_rational, 2> const& polynomial) {
        return function_scale * Compose(polynomial,
                                        AccuratePolynomial<cpp_rational, 1>(
                                            {0, 1 / argument_scale}));
      };
  return {.functions = scaled_functions,
          .polynomials = {build_scaled_polynomial(function_scales[0],
                                                  polynomials[0]),
                          build_scaled_polynomial(function_scales[1],
                                                  polynomials[1])},
          .remainders = scaled_remainders,
          .argument = scaled_argument};
}

// This is essentially the same as Gal's exhaustive search, but with the
// scaling to binade 0.
absl::StatusOr<std::int64_t> StehléZimmermannExhaustiveSearch(
    std::array<AccurateFunction, 2> const& F,
    std::int64_t const M,
    std::int64_t const T) {
  VLOG(3) << "Exhaustive search with T = " << T;
  for (std::int64_t t = 0; t <= T; ++t) {
    {
      bool found = true;
      for (auto const& Fᵢ : F) {
        auto const Fᵢ_t = Fᵢ(t);
        auto const Fᵢ_t_cmod_1 = Fᵢ_t - round(Fᵢ_t);
        VLOG(3) << "Fi(t) cmod 1 = " << Fᵢ_t_cmod_1;
        if (M * abs(Fᵢ_t_cmod_1) >= 1) {
          found = false;
          break;
        }
      }
      if (found) {
        VLOG(3) << "t = " << t;
        return t;
      }
    }
    if (t > 0) {
      bool found = true;
      for (auto const& Fᵢ : F) {
        auto const Fᵢ_minus_t = Fᵢ(-t);
        auto const Fᵢ_minus_t_cmod_1 = Fᵢ_minus_t - round(Fᵢ_minus_t);
        VLOG(3) << "Fi(-t) cmod 1 = " << Fᵢ_minus_t_cmod_1;
        if (M * abs(Fᵢ_minus_t_cmod_1) >= 1) {
          found = false;
          break;
        }
      }
      if (found) {
        VLOG(3) << "t = " << -t;
        return -t;
      }
    }
  }
  return absl::NotFoundError("Not enough zeroes");
}

// Searches in a "slice", which is a set of two intervals of measure `2 * T₀` on
// either side of `scaled.argument`.  Consecutive values of `slice_index`
// correspond to contiguous intervals farther away from `scaled.argument`.
// Slices may be processed independently of one another.
// Returns a *scaled* argument, or `NotFound` if no solution was found in the
// slice.
template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousSliceSearch(
    StehléZimmermannSpecification const& scaled,
    std::int64_t const slice_index) {
  constexpr std::int64_t M = 1LL << zeroes;
  constexpr std::int64_t N = 1LL << std::numeric_limits<double>::digits;

  // [SZ05], section 3.2, proves that T³ = O(M * N).  Of course, the
  // multiplicative factor is not known.  In practice it seems that a value of
  // T₀ that's too large is very costly as it results in many intervals that are
  // rejected with `OutOfRange` and must be halved and retried.  A value that's
  // too small on the other hand can slow down progress.  The fudge factor 1/128
  // attempts to strike a balance between these problems; it has been chosen by
  // benchmarking SinCos18 around 1167/2048.
  std::int64_t const T₀ = static_cast<std::int64_t>(
      Cbrt(static_cast<double>(M) * static_cast<double>(N)) / 128.0);

  // Construct intervals of measure `2 * T₀` above and below `scaled.argument`
  // and search for solutions on each side alternatively.
  Interval<cpp_rational> const initial_high_interval{
      .min = scaled.argument + cpp_rational(2 * slice_index * T₀, N),
      .max = scaled.argument + cpp_rational(2 * (slice_index + 1) * T₀, N)};
  Interval<cpp_rational> const initial_low_interval{
      .min = scaled.argument - cpp_rational(2 * (slice_index + 1) * T₀, N),
      .max = scaled.argument - cpp_rational(2 * slice_index * T₀, N)};

  Interval<cpp_rational> high_interval = initial_high_interval;
  Interval<cpp_rational> low_interval = initial_low_interval;

  // The radii of the intervals remaining to cover above and below the
  // `scaled.argument`.
  std::int64_t high_T_to_cover = T₀;
  std::int64_t low_T_to_cover = T₀;

  // When exiting this loop, we have completely processed
  // `initial_high_interval` and `initial_low_interval`.
  for (;;) {
    bool const high_interval_empty = high_interval.empty();
    bool const low_interval_empty = low_interval.empty();
    if (high_interval_empty && low_interval_empty) {
      return absl::NotFoundError(
          absl::StrCat("No solution in slice #", slice_index));
    }

    if (!high_interval_empty) {
      std::int64_t T = high_T_to_cover;
      // This loop exits (breaks or returns) when `T <= T_max` because
      // exhaustive search always gives an answer.
      for (;;) {
        VLOG(3) << "T = " << T << ", high_interval = " << high_interval;
        auto const status_or_solution =
            StehléZimmermannSimultaneousSearch<zeroes>(scaled.functions,
                                                       scaled.polynomials,
                                                       scaled.remainders,
                                                       high_interval.midpoint(),
                                                       N,
                                                       T);
        absl::Status const& status = status_or_solution.status();
        if (status.ok()) {
          return status_or_solution.value();
        } else {
          VLOG(3) << "Status = " << status;
          if (absl::IsOutOfRange(status)) {
            // Halve the interval.  Make sure that the new interval is
            // contiguous to the segment already explored.
            T /= 2;
            high_interval.max = high_interval.min + cpp_rational(2 * T, N);
          } else if (absl::IsNotFound(status)) {
            // No solutions here, go to the next interval.
            high_T_to_cover -= T;
            break;
          } else {
            return status;
          }
        }
      }
    }
    if (!low_interval_empty) {
      std::int64_t T = low_T_to_cover;
      // This loop exits (breaks or returns) when `T <= T_max` because
      // exhaustive search always gives an answer.
      for (;;) {
        VLOG(3) << "T = " << T << ", low_interval = " << low_interval;
        auto const status_or_solution =
            StehléZimmermannSimultaneousSearch<zeroes>(scaled.functions,
                                                       scaled.polynomials,
                                                       scaled.remainders,
                                                       low_interval.midpoint(),
                                                       N,
                                                       T);
        absl::Status const& status = status_or_solution.status();
        if (status.ok()) {
          return status_or_solution.value();
        } else {
          VLOG(3) << "Status = " << status;
          if (absl::IsOutOfRange(status)) {
            // Halve the interval.  Make sure that the new interval is
            // contiguous to the segment already explored.
            T /= 2;
            low_interval.min = low_interval.max - cpp_rational(2 * T, N);
          } else if (absl::IsNotFound(status)) {
            // No solutions here, go to the next interval.
            low_T_to_cover -= T;
            break;
          } else {
            return status;
          }
        }
      }
    }
    VLOG_EVERY_N(2, 10) << "high = "
                        << DebugString(static_cast<double>(high_interval.max));
    VLOG_EVERY_N(2, 10) << "low  = "
                        << DebugString(static_cast<double>(low_interval.min));
    high_interval = {.min = high_interval.max,
                     .max = initial_high_interval.max};
    low_interval = {.min = initial_low_interval.min, .max = low_interval.min};
  }
}


template<std::int64_t zeroes>
cpp_rational GalExhaustiveSearch(std::vector<AccurateFunction> const& functions,
                                 cpp_rational const& starting_argument) {
  CHECK_LT(0, starting_argument);

  // We will look for candidates both above and below `starting_argument`.
  // Note that if `starting_argument` is a power of 2, the increments above
  // and below `starting_argument` are not the same.
  std::int64_t exponent;
  auto const starting_mantissa =
      frexp(static_cast<cpp_bin_float_50>(starting_argument), &exponent);
  cpp_rational const high_increment =
      exp2(exponent - std::numeric_limits<double>::digits);
  cpp_rational const low_increment =
      starting_mantissa == 0.5 ? high_increment / 2 : high_increment;

  cpp_rational high_x = starting_argument;
  cpp_rational low_x = starting_argument - low_increment;
  for (;;) {
    if (AllFunctionValuesHaveDesiredZeroes<zeroes>(functions, high_x)) {
      return high_x;
    }
    high_x += high_increment;
    if (AllFunctionValuesHaveDesiredZeroes<zeroes>(functions, low_x)) {
      return low_x;
    }
    low_x -= low_increment;
  }
}

template<std::int64_t zeroes>
std::vector<cpp_rational> GalExhaustiveMultisearch(
    std::vector<AccurateFunction> const& functions,
    std::vector<cpp_rational> const& starting_arguments) {
  ThreadPool<cpp_rational> search_pool(std::thread::hardware_concurrency());

  std::vector<std::future<cpp_rational>> futures;
  for (std::int64_t i = 0; i < starting_arguments.size(); ++i) {
    futures.push_back(search_pool.Add([i, &functions, &starting_arguments]() {
      return GalExhaustiveSearch<zeroes>(functions, starting_arguments[i]);
    }));
  }

  std::vector<cpp_rational> results;
  for (auto& future : futures) {
    results.push_back(future.get());
  }
  return results;
}

template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousSearch(
    std::array<AccurateFunction, 2> const& functions,
    std::array<AccuratePolynomial<cpp_rational, 2>, 2> const& polynomials,
    std::array<AccurateFunction, 2> const& remainders,
    cpp_rational const& starting_argument,
    std::int64_t const N,
    std::int64_t const T) {
  // This implementation follows [SZ05], section 3.1.
  std::int64_t const M = 1 << zeroes;

  // Preliminary: shift and rescale the functions and the polynomials:
  //   Fᵢ = N fᵢ(t / N)
  //   Pᵢ = N pᵢ(t / N)
  std::array<AccurateFunction, 2> F;
  std::array<std::optional<AccuratePolynomial<cpp_rational, 2>>, 2> P;
  for (std::int64_t i = 0; i < functions.size(); ++i) {
    F[i] = [&functions, i, N, &starting_argument](cpp_rational const& t) {
      // Here |t| ≤ T.
      return N * functions[i](starting_argument + t / N);
    };
  }

  // If the interval is small enough, we don't use the Stehlé-Zimmermann
  // algorithm.  Instead, we use an exhaustive search.  Note that this may yield
  // a better solution, because if there is one in the interval, it is sure to
  // find it, whereas Stehlé-Zimmermann may miss it.
  if (T <= T_max) {
    auto const status_or_t = StehléZimmermannExhaustiveSearch(F, M, T);
    RETURN_IF_ERROR(status_or_t);

    return starting_argument + cpp_rational(status_or_t.value(), N);
  }

  AccuratePolynomial<cpp_rational, 1> const shift_and_rescale(
      {starting_argument, cpp_rational(1, N)});
  for (std::int64_t i = 0; i < polynomials.size(); ++i) {
    P[i] = N * Compose(polynomials[i], shift_and_rescale);
  }

  // Step 2: compute ε.  We use the remainders provided by the clients.  Note
  // that we could save on the number of evaluations by providing both bounds to
  // a single call.
  cpp_bin_float_50 ε = 0;
  for (std::int64_t i = 0; i < remainders.size(); ++i) {
    auto const T_over_N = cpp_rational(T, N);
    ε = std::max(ε, abs(N * remainders[i](starting_argument - T_over_N)));
    ε = std::max(ε, abs(N * remainders[i](starting_argument + T_over_N)));
  }
  VLOG(3) << "ε = " << ε;

  // Step 3, first part: compute Mʹ and C.  Give up is C is 0, which may happen
  // if ε is too large.
  auto const Mʹ = static_cast<std::int64_t>(floor(M / (2 + 2 * M * ε)));
  auto const C = 3 * Mʹ;
  if (C == 0) {
    return absl::FailedPreconditionError("Error too large");
  }
  VLOG(3) << "C = " << C;

  // Step 3, second part: compute P̃
  std::array<std::optional<AccuratePolynomial<cpp_int, 2>>, 2> P̃;
  AccuratePolynomial<cpp_rational, 1> const Tτ({0, T});
  for (std::int64_t i = 0; i < P.size(); ++i) {
    auto const composition_coefficients = Compose(C * *P[i], Tτ).coefficients();
    AccuratePolynomial<cpp_int, 2>::Coefficients P̃_coefficients;
    for_all_of(composition_coefficients, P̃_coefficients)
        .loop([](auto const& composition_coefficient, auto& P̃_coefficient) {
          P̃_coefficient = static_cast<cpp_int>(Round(composition_coefficient));
        });
    P̃[i] = AccuratePolynomial<cpp_int, 2>(P̃_coefficients);
    VLOG(3) << "P̃[" << i << "] = " << *P̃[i];
  }

  // Step 5 and 6: form the lattice.  Note that our vectors are in columns, not
  // in rows as in the paper.
  auto const& P̃₀_coefficients = P̃[0]->coefficients();
  auto const& P̃₁_coefficients = P̃[1]->coefficients();

  using Lattice = FixedMatrix<cpp_int, 5, 4>;

  Lattice const L(
      {C,     0, std::get<0>(P̃₀_coefficients), std::get<0>(P̃₁_coefficients),
       0, C * T, std::get<1>(P̃₀_coefficients), std::get<1>(P̃₁_coefficients),
       0,     0, std::get<2>(P̃₀_coefficients), std::get<2>(P̃₁_coefficients),
       0,     0,                            3,                            0,
       0,     0,                            0,                            3});
  VLOG(3) << "L = " << L;

  // Step 7: reduce the lattice.
  // The lattice really has integer coefficients, but this is inconvenient to
  // propagate through the matrix algorithms.  (It would require copies instead
  // of views for all the types, not just the ones we use here.)
  Lattice const V = NguyễnStehlé(L);
  VLOG(3) << "V = " << V;

  // Step 8: find the three shortest vectors of the reduced lattice.  We sort
  // the columns according to the L₂ norm.
  std::array<std::unique_ptr<ColumnView<Lattice const>>, V.columns()> v;
  for (std::int64_t i = 0; i < v.size(); ++i) {
    v[i] = std::make_unique<ColumnView<Lattice const>>(
        ColumnView<Lattice const>{.matrix = V,
                                  .first_row = 0,
                                  .last_row = V.rows() - 1,
                                  .column = i});
  }
  std::sort(v.begin(),
            v.end(),
            [](std::unique_ptr<ColumnView<Lattice const>> const& left,
               std::unique_ptr<ColumnView<Lattice const>> const& right) {
              return left->Norm²() < right->Norm²();
            });

  // Step 9: check that the shortest vectors are short enough.
  auto norm1 = [](ColumnView<Lattice const> const& v) {
    cpp_int norm1 = 0;
    for (std::int64_t i = 0; i < v.size(); ++i) {
      norm1 += abs(v[i]);
    }
    return norm1;
  };

  static constexpr std::int64_t dimension = 3;
  for (std::int64_t i = 0; i < dimension; ++i) {
    auto const& vᵢ = *v[i];
    VLOG(3) << "v[" << i << "] = " << vᵢ;
    if (norm1(vᵢ) >= C) {
      return absl::OutOfRangeError("Vectors too big");
    }
  }

  // Step 10: compute Q by eliminating the last two variables.  Give up if the
  // degree 1 coefficient is 0, there is no solution.
  std::array<cpp_int, dimension> Q_multipliers;
  for (std::int64_t i = 0; i < dimension; ++i) {
    auto const& vᵢ₊₁ = *v[(i + 1) % dimension];
    auto const& vᵢ₊₂ = *v[(i + 2) % dimension];
    Q_multipliers[i] = vᵢ₊₁[3] * vᵢ₊₂[4] - vᵢ₊₁[4] * vᵢ₊₂[3];
  }

  FixedVector<cpp_int, 2> Q_coefficients{};
  for (std::int64_t i = 0; i < dimension; ++i) {
    auto const& vᵢ = *v[i];
    for (std::int64_t j = 0; j < Q_coefficients.size(); ++j) {
      Q_coefficients[j] += Q_multipliers[i] * vᵢ[j];
    }
  }

  AccuratePolynomial<cpp_rational, 1> const Q({Q_coefficients[0],
                                               Q_coefficients[1]});
  VLOG(3) << "Q = " << Q;
  if (Q_coefficients[1] == 0) {
      return absl::NotFoundError("No integer zeroes");
  }

  // Step 11: compute q and find its integer root (singular), if any.
  AccuratePolynomial<cpp_rational, 1> const q =
      Compose(Q, AccuratePolynomial<cpp_rational, 1>({0, cpp_rational(1, T)}));

  cpp_rational const t₀ =
      -std::get<0>(q.coefficients()) / std::get<1>(q.coefficients());
  VLOG(3) << "t₀ = " << t₀;
  if (abs(t₀) > T) {
    return absl::NotFoundError("Out of bounds");
  } else if (denominator(t₀) != 1) {
    return absl::NotFoundError("Noninteger root");
  }

  for (auto const& Fᵢ : F) {
    auto const Fᵢ_t₀ = Fᵢ(t₀);
    auto const Fᵢ_t₀_cmod_1 = Fᵢ_t₀ - round(Fᵢ_t₀);
    VLOG(3) << "Fi(t₀) cmod 1 = " << Fᵢ_t₀_cmod_1;
    if (M * abs(Fᵢ_t₀_cmod_1) >= 1) {
      return absl::NotFoundError("Not enough zeroes");
    }
  }

  return starting_argument + t₀ / N;
}

template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousFullSearch(
    std::array<AccurateFunction, 2> const& functions,
    std::array<AccuratePolynomial<cpp_rational, 2>, 2> const& polynomials,
    std::array<AccurateFunction, 2> const& remainders,
    cpp_rational const& starting_argument,
    ThreadPool<void>* const search_pool) {
  // Start by scaling the specification of the search.  The rest of this
  // function only uses the scaled objects.
  double argument_scale;
  auto const scaled = ScaleToBinade0({.functions = functions,
                                      .polynomials = polynomials,
                                      .remainders = remainders,
                                      .argument = starting_argument},
                                     argument_scale);

  // This mutex is not contended as it is only held exclusively when we have
  // found a solution.
  absl::Mutex lock;
  std::optional<absl::StatusOr<cpp_rational>> status_or_solution
      GUARDED_BY(lock);

  auto search_one_slice = [argument_scale,
                           &lock,
                           &scaled,
                           &starting_argument,
                           &status_or_solution](
                              std::int64_t const slice_index) {
    auto const start = std::chrono::system_clock::now();

    auto const status_or_scaled_solution =
        StehléZimmermannSimultaneousSliceSearch<zeroes>(scaled, slice_index);
    auto const end = std::chrono::system_clock::now();
    VLOG(1) << "Search for slice #" << slice_index << " around "
            << starting_argument << " took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                     start);

    absl::Status const& status = status_or_scaled_solution.status();
    if (status.ok()) {
      // The argument returned by the slice search is scaled, so we must
      // adjust it before returning.
      auto const solution = status_or_scaled_solution.value() / argument_scale;
      VLOG(1) << "Solution for " << starting_argument << ", slice #"
              << slice_index;
      {
        absl::MutexLock l(&lock);
        // We have found a solution; we only retain it if (1) no internal error
        // occurred; and (2) it closer to the `starting_argument` than any
        // solution found previously.
        if (status_or_solution.has_value()) {
          if (status_or_solution.value().ok() &&
              abs(solution - starting_argument) <
                  abs(status_or_solution.value().value() - starting_argument)) {
            status_or_solution = solution;
          }
        } else {
          status_or_solution = solution;
        }
      }
    } else if (absl::IsNotFound(status)) {
      // No solution found in this slice, go to the next one.
    } else {
      // Some kind of internal error, give up.
      {
        absl::MutexLock l(&lock);
        status_or_solution = status;
      }
    }
  };

  // This variable should not be under a mutex as it might cause contention.
  // Apparently the memory barrier implied by sequential consistency is not
  // degrading performance.
  std::atomic<std::int64_t> current_slice_index = 0;

  // This thread attempts to keep CPU utilization at 100% by starting
  // speculative searches if some of the threads in the `search_pool` are idle.
  // When there are calls queued in the `search_pool`, it does essentially
  // nothing, idly looping every 1 s.
  std::vector<std::future<void>> speculative_futures;
  std::thread speculative_scheduler([&current_slice_index,
                                     &lock,
                                     &search_one_slice,
                                     search_pool,
                                     &speculative_futures,
                                     &starting_argument,
                                     &status_or_solution]() {
    if (search_pool != nullptr) {
      for (;;) {
        // Avoid busy waiting, and only wake up if the pool has an idle thread,
        // or if we have been stuck for too long.  The timeout ensures that this
        // loop eventually terminates if a solution is found without speculative
        // execution.  Of course, this is racy, so there is no guarantee that
        // the next call to `TryAdd` below will succeed.
        search_pool->WaitUntilIdleFor(absl::Seconds(1));

        // Stop this thread if a solution has been found.
        {
          absl::ReaderMutexLock l(&lock);
          if (status_or_solution.has_value()) {
            return;
          }
        }

        // Try to queue as many speculative searches as possible to keep the
        // `search_pool` busy.  Note that the slice to work on is determined
        // when the function actually starts.
        for (;;) {
          auto maybe_future = search_pool->TryAdd(
              [&search_one_slice, &current_slice_index, &starting_argument] {
                std::int64_t const slice_index =
                    current_slice_index.fetch_add(1);
                VLOG(1) << "Speculative search for " << starting_argument
                        << ", slice #" << slice_index;
                search_one_slice(slice_index);
              });
          if (!maybe_future.has_value()) {
            break;
          }
          speculative_futures.push_back(std::move(maybe_future).value());
        }
      }
    }
  });

  for (;;) {
    VLOG(1) << "Sequential search for " << starting_argument << ", slice #"
            << current_slice_index;
    search_one_slice(current_slice_index.fetch_add(1));

    absl::ReaderMutexLock l(&lock);
    if (status_or_solution.has_value()) {
      break;
    }
  }

  // Wait for any remaining speculative execution to complete.  They may find a
  // better solution than the one that caused us to exit the sequential loop.
  for (auto const& future : speculative_futures) {
    future.wait();
  }
  speculative_scheduler.join();

  return status_or_solution.value();
}

template<std::int64_t zeroes>
std::vector<absl::StatusOr<cpp_rational>>
StehléZimmermannSimultaneousMultisearch(
    std::array<AccurateFunction, 2> const& functions,
    std::vector<std::array<AccuratePolynomial<cpp_rational, 2>, 2>> const&
        polynomials,
    std::vector<std::array<AccurateFunction, 2>> const& remainders,
    std::vector<cpp_rational> const& starting_arguments) {
  std::vector<absl::StatusOr<cpp_rational>> result;
  result.resize(starting_arguments.size());
  StehléZimmermannSimultaneousStreamingMultisearch<zeroes>(
      functions,
      polynomials,
      remainders,
      starting_arguments,
      [&result](std::int64_t const index,
                absl::StatusOr<cpp_rational> status_or_final_argument) {
        result[index] = std::move(status_or_final_argument);
      });
  return result;
}

template<std::int64_t zeroes>
void StehléZimmermannSimultaneousStreamingMultisearch(
    std::array<AccurateFunction, 2> const& functions,
    std::vector<std::array<AccuratePolynomial<cpp_rational, 2>, 2>> const&
        polynomials,
    std::vector<std::array<AccurateFunction, 2>> const& remainders,
    std::vector<cpp_rational> const& starting_arguments,
    std::function<void(/*index=*/std::int64_t,
                       absl::StatusOr<cpp_rational>)> const& callback) {
  ThreadPool<void> search_pool(std::thread::hardware_concurrency());

  std::vector<std::future<void>> futures;
  for (std::int64_t i = 0; i < starting_arguments.size(); ++i) {
    futures.push_back(search_pool.Add([i,
                                       &callback,
                                       &functions,
                                       &polynomials,
                                       &remainders,
                                       &search_pool,
                                       &starting_arguments]() {
      auto const& starting_argument = starting_arguments[i];
      LOG(INFO) << "Starting search around " << starting_argument;
      auto status_or_final_argument =
          StehléZimmermannSimultaneousFullSearch<zeroes>(functions,
                                                         polynomials[i],
                                                         remainders[i],
                                                         starting_argument,
                                                         &search_pool);
      if (status_or_final_argument.ok()) {
        LOG(INFO) << "Finished search around " << starting_argument
                  << ", found " << status_or_final_argument.value();
      } else {
        LOG(WARNING) << "Search around " << starting_argument << " failed with"
                     << status_or_final_argument.status();
      }
      callback(i, std::move(status_or_final_argument));
    }));
  }

  for (auto const& future : futures) {
    future.wait();
  }
}

}  // namespace internal
}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
