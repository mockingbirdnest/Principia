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

// For intervals with a radius less or equal to this value, we use exhaustive
// search.
constexpr std::int64_t T_max = 4;
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
  std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2> polynomials;
  std::array<ApproximateFunctionFactory, 2> remainders;
  cpp_rational argument;
};

template<typename Factory>
std::array<std::invoke_result_t<Factory, cpp_rational>, 2> EvaluateFactoriesAt(
    std::array<Factory, 2> const& factories,
    cpp_rational const& argument) {
  return {factories[0](argument), factories[1](argument)};
}

// In general, scales the argument, functions, polynomials, and remainders to
// lie within [1/2, 1[.  There is a subtlety though if the input is such that
// either the argument or a function is a power of 2 or close enough to a power
// of 2 that it could cross from one binade to another as part of the search.
// The Stehlé-Zimmermann algorithm wants to believe that machine numbers are
// equally spaced, which is not the case when changing binade.  We have two
// options:
// 1. Choose the scale based on the spacing near 0, which results in a larger
//    scale factor and can possibly cause bogus solutions to be found near ∞
//    (such solutions would be midpoints between machine numbers and not machine
//    numbers).  This may cause the argument and functions to lie within
//    [1/2, 2[.
// 2. Choose the scale based on the spacing near ∞, which results in a smaller
//    scale factor and can possibly cause solutions to be missed near 0 because
//    the algorithm would skip half of the machine numbers.
// Option 1 is the only one that ensures optimality of the solution, but it
// requires that each candidate solution be checked to see if it actually
// fullfills the desired conditions.
template<std::int64_t zeroes>
StehléZimmermannSpecification ScaleToBinade01(
    StehléZimmermannSpecification const& input,
    double& argument_scale) {
  auto const& functions = input.functions;
  auto const& polynomials = input.polynomials;
  auto const& remainders = input.remainders;
  auto const& starting_argument = input.argument;

  // In order to decide whether there is a risk of changing binade, we need an
  // estimate of the interval that will be searched.  We use the random model of
  // [SZ05], section 1, together with a safety margin in case the search would
  // be unlucky.
  constexpr std::int64_t multiplier_safety_bits = 6;
  constexpr double multiplier =
      1 + std::numeric_limits<double>::epsilon() *
              (1LL << (2 * zeroes - 2 + multiplier_safety_bits));
  auto const lower_bound = starting_argument / multiplier;
  auto const upper_bound = starting_argument * multiplier;

  // Returns a scale factor suitable for both `value1` and `value2`.  Makes no
  // assumptions regarding the order of its arguments.
  auto const compute_scale = [](auto const& value1,
                                auto const& value2) {
    std::int64_t value1_exponent;
    auto const value1_mantissa =
        frexp(static_cast<cpp_bin_float_50>(value1), &value1_exponent);
    CHECK_NE(0, value1_mantissa);

    std::int64_t value2_exponent;
    auto const value2_mantissa =
        frexp(static_cast<cpp_bin_float_50>(value2), &value2_exponent);
    CHECK_NE(0, value2_mantissa);

    CHECK_LE(std::abs(value1_exponent - value2_exponent), 1)
        << "Values straddle multiple powers of 2: " << value1 << " and "
        << value2;

    return exp2(std::max(-value1_exponent, -value2_exponent));
  };

  argument_scale = compute_scale(lower_bound, upper_bound);
  cpp_rational const scaled_argument = starting_argument * argument_scale;

  std::array<double, 2> function_scales;
  std::array<AccurateFunction, 2> scaled_functions;
  std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2> scaled_polynomials;
  std::array<ApproximateFunctionFactory, 2> scaled_remainders;
  for (std::int64_t i = 0; i < scaled_functions.size(); ++i) {
    function_scales[i] = compute_scale(functions[i](lower_bound),
                                       functions[i](upper_bound));
    scaled_functions[i] = [argument_scale,
                           function_scale = function_scales[i],
                           function =
                               functions[i]](cpp_rational const& argument) {
      return function_scale * function(argument / argument_scale);
    };
    scaled_polynomials[i] = [argument_scale,
                             function_scale = function_scales[i],
                             polynomial = polynomials[i]](
                                cpp_rational const& argument₀) {
      return function_scale * Compose(polynomial(argument₀ / argument_scale),
                                      AccuratePolynomial<cpp_rational, 1>(
                                          {0, 1 / argument_scale}));
    };
    scaled_remainders[i] = [argument_scale,
                            function_scale = function_scales[i],
                            remainder =
                                remainders[i]](cpp_rational const& argument₀) {
      return [argument₀,
              argument_scale,
              function_scale,
              remainder](cpp_rational const& argument) {
        return function_scale *
               remainder(argument₀ / argument_scale)(argument / argument_scale);
      };
    };
  }

  return {.functions = scaled_functions,
          .polynomials = scaled_polynomials,
          .remainders = scaled_remainders,
          .argument = scaled_argument};
}

// `ScaleToBinade01` may have had to choose a scale factor that's "too large" to
// account for the argument or functions possibly changing binade.  This may
// result in candidate solutions that are midway between machine numbers.  This
// function rejects such solutions.
template<std::int64_t zeroes>
bool VerifyBinade01Solution(StehléZimmermannSpecification const& scaled,
                            cpp_rational const& scaled_solution) {
  constexpr std::int64_t M = 1LL << zeroes;
  CHECK_LE(cpp_rational(1, 2), abs(scaled_solution));

  // If the scaled solution is below 1, it is necessarily a machine number.  It
  // may be above 1 if we upscaled, in which case we must verify that once
  // scaled down to binade 0 it is a machine number.
  if (abs(scaled_solution) > 1) {
    auto const solution_binade0 = scaled_solution / 2;
    if (solution_binade0 != static_cast<double>(solution_binade0)) {
      VLOG(1) << "Rejecting " << scaled_solution
              << " because it's not a machine number";
      return false;
    }
  }

  // Verify that the value of the functions are close to a machine number.  This
  // may not be the case if we upscaled the functions.
  for (auto const& function : scaled.functions) {
    auto const y = function(scaled_solution);
    std::int64_t y_exponent;
    auto const y_mantissa =
        frexp(static_cast<cpp_bin_float_50>(y), &y_exponent);
    auto const y_mantissa_scaled =
        ldexp(y_mantissa, std::numeric_limits<double>::digits);
    auto const y_mantissa_scaled_cmod_1 =
        y_mantissa_scaled - round(y_mantissa_scaled);
    if (M * abs(y_mantissa_scaled_cmod_1) >= 1) {
      VLOG(1) << "Rejecting " << scaled_solution << " because " << y
              << " is not close enough to a machine number";
      return false;
    }
  }

  return true;
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
  // multiplicative factor is not known, but 1 works well in practice.
  std::int64_t const T₀ = static_cast<std::int64_t>(
      Cbrt(static_cast<double>(M) * static_cast<double>(N)));

  // Construct intervals of measure `2 * T₀` above and below `scaled.argument`
  // and search for solutions on each side alternatively.
  Interval<cpp_rational> const initial_high_interval{
      .min = scaled.argument + cpp_rational(2 * slice_index * T₀, N),
      .max = scaled.argument + cpp_rational(2 * (slice_index + 1) * T₀, N)};
  Interval<cpp_rational> const initial_low_interval{
      .min = scaled.argument - cpp_rational(2 * (slice_index + 1) * T₀, N),
      .max = scaled.argument - cpp_rational(2 * slice_index * T₀, N)};

  // Evaluate the factories at the centre of each half of the slice.
  auto const high_polynomials =
      EvaluateFactoriesAt(scaled.polynomials, initial_high_interval.midpoint());
  auto const high_remainders =
      EvaluateFactoriesAt(scaled.remainders, initial_high_interval.midpoint());
  auto const low_polynomials =
      EvaluateFactoriesAt(scaled.polynomials, initial_low_interval.midpoint());
  auto const low_remainders =
      EvaluateFactoriesAt(scaled.remainders, initial_low_interval.midpoint());

  // The radii of the intervals remaining to cover above and below the
  // `scaled.argument`.
  std::int64_t high_T_to_cover = T₀;
  std::int64_t low_T_to_cover = T₀;

  // When exiting this loop, we have completely processed
  // `initial_high_interval` and `initial_low_interval`.
  for (;;) {
    if (high_T_to_cover == 0 && low_T_to_cover == 0) {
      return absl::NotFoundError(
          absl::StrCat("No solution in slice #", slice_index));
    }

    if (high_T_to_cover > 0) {
      std::int64_t T = high_T_to_cover;
      // This loop exits (breaks or returns) when `T <= T_max` because
      // exhaustive search always gives an answer.
      for (;;) {
        // Make sure that the new interval is contiguous to the segment already
        // explored.
        cpp_rational const high_interval_midpoint =
            initial_high_interval.max -
            cpp_rational(2 * high_T_to_cover - T, N);
        VLOG(3) << "T = " << T << ", high_T_to_cover = " << high_T_to_cover;
        auto const status_or_solution =
            StehléZimmermannSimultaneousSearch<zeroes>(scaled.functions,
                                                       high_polynomials,
                                                       high_remainders,
                                                       high_interval_midpoint,
                                                       N,
                                                       T);
        absl::Status const& status = status_or_solution.status();
        if (status.ok()) {
          auto const& solution = status_or_solution.value();
          if (VerifyBinade01Solution<zeroes>(scaled, solution)) {
            return solution;
          } else if (T == 1) {
            high_T_to_cover -= T;
          } else {
            T /= 2;
          }
        } else {
          VLOG(3) << "Status = " << status;
          if (absl::IsOutOfRange(status)) {
            // Halve the interval.
            T /= 2;
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
    if (low_T_to_cover > 0) {
      std::int64_t T = low_T_to_cover;
      // This loop exits (breaks or returns) when `T <= T_max` because
      // exhaustive search always gives an answer.
      for (;;) {
        // Make sure that the new interval is contiguous to the segment already
        // explored.
        cpp_rational const low_interval_midpoint =
            initial_low_interval.min +
            cpp_rational(2 * low_T_to_cover - T, N);
        VLOG(3) << "T = " << T << ", low_T_to_cover = " << low_T_to_cover;
        auto const status_or_solution =
            StehléZimmermannSimultaneousSearch<zeroes>(scaled.functions,
                                                       low_polynomials,
                                                       low_remainders,
                                                       low_interval_midpoint,
                                                       N,
                                                       T);
        absl::Status const& status = status_or_solution.status();
        if (status.ok()) {
          auto const& solution = status_or_solution.value();
          if (VerifyBinade01Solution<zeroes>(scaled, solution)) {
            return solution;
          } else if (T == 1) {
            low_T_to_cover -= T;
          } else {
            T /= 2;
          }
        } else {
          VLOG(3) << "Status = " << status;
          if (absl::IsOutOfRange(status)) {
            // Halve the interval.  Make sure that the new interval is
            // contiguous to the segment already explored.
            T /= 2;
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
    std::array<ApproximateFunction, 2> const& remainders,
    cpp_rational const& starting_argument,
    std::int64_t const N,
    std::int64_t const T) {
  // This implementation follows [SZ05], section 3.1.
  std::int64_t const M = 1 << zeroes;

  // Preliminary: shift and rescale the functions and the polynomials:
  //   Fᵢ = N fᵢ(t / N)
  //   Pᵢ = N pᵢ(t / N)
  // We don't reify Pᵢ as that would require an extra call to `Compose`.
  // Instead, we directly compute P̃ᵢ below.
  std::array<AccurateFunction, 2> F;
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

  // Step 2: compute ε.  We use the remainders provided by the clients.  Note
  // that we could save on the number of evaluations by providing both bounds to
  // a single call.
  double ε = 0;
  auto const T_over_N = cpp_rational(T, N);
  for (std::int64_t i = 0; i < remainders.size(); ++i) {
    ε = std::max(ε, std::abs(N * remainders[i](starting_argument - T_over_N)));
    ε = std::max(ε, std::abs(N * remainders[i](starting_argument + T_over_N)));
  }
  VLOG(3) << "ε = " << ε;

  // Step 3, first part: compute Mʹ and C.  Give up is C is 0, which may happen
  // if ε is too large.
  auto const Mʹ = static_cast<std::int64_t>(std::floor(M / (2 + 2 * M * ε)));
  auto const C = 3 * Mʹ;
  if (C == 0) {
    return absl::FailedPreconditionError("Error too large");
  }
  VLOG(3) << "C = " << C;

  // Step 3, second part: compute P̃
  std::array<std::optional<AccuratePolynomial<cpp_int, 2>>, 2> P̃;
  AccuratePolynomial<cpp_rational, 1> const shift_and_rescale(
      {starting_argument, cpp_rational(T, N)});
  for (std::int64_t i = 0; i < polynomials.size(); ++i) {
    auto const composition_coefficients =
        Compose(C * (N * polynomials[i]), shift_and_rescale).coefficients();
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
      LOG_IF(FATAL, Q_coefficients[0] == 0) << "Identically zero";
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
    std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2> const&
        polynomials,
    std::array<ApproximateFunctionFactory, 2> const& remainders,
    cpp_rational const& starting_argument,
    ThreadPool<void>* const search_pool) {
  // Start by scaling the specification of the search.  The rest of this
  // function only uses the scaled objects.
  double argument_scale;
  auto const scaled = ScaleToBinade01<zeroes>({.functions = functions,
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
      absl::MutexLock l(&lock);
      // The argument returned by the slice search is scaled, so we must
      // adjust it before returning.
      auto const scaled_solution = status_or_scaled_solution.value();
      auto const solution = scaled_solution / argument_scale;
      VLOG(1) << "Solution for " << starting_argument << ", slice #"
              << slice_index << " is " << solution;
      // We have found a solution; we only retain it if (1) no internal error
      // occurred; and (2) it closer to the `starting_argument` than any
      // solution found previously.
      if (status_or_solution.has_value()) {
        if (status_or_solution.value().ok()) {
          if (abs(solution - starting_argument) <
              abs(status_or_solution.value().value() - starting_argument)) {
            status_or_solution = solution;
          } else {
            VLOG(1) << "Solution for slice #" << slice_index
                    << " discarded because there is a better one";
          }
        }
      } else {
        status_or_solution = solution;
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
    std::int64_t const slice_index = current_slice_index.fetch_add(1);
    VLOG(1) << "Sequential search for " << starting_argument << ", slice #"
            << slice_index;
    search_one_slice(slice_index);

    absl::ReaderMutexLock l(&lock);
    if (status_or_solution.has_value()) {
      break;
    }
  }

  // Wait for any remaining speculative execution to complete.  They may find a
  // better solution than the one that caused us to exit the sequential loop.
  // Note that it's important to join the scheduler first, so that no more
  // speculative work is started.
  speculative_scheduler.join();
  for (auto const& future : speculative_futures) {
    future.wait();
  }

  return status_or_solution.value();
}

template<std::int64_t zeroes>
std::vector<absl::StatusOr<cpp_rational>>
StehléZimmermannSimultaneousMultisearch(
    std::array<AccurateFunction, 2> const& functions,
    std::vector<std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2>>
        const& polynomials,
    std::vector<std::array<ApproximateFunctionFactory, 2>> const& remainders,
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
    std::vector<std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2>>
        const& polynomials,
    std::vector<std::array<ApproximateFunctionFactory, 2>> const& remainders,
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
