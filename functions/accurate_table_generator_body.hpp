#pragma once

#include "functions/accurate_table_generator.hpp"

#include <algorithm>
#include <future>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/tags.hpp"
#include "base/thread_pool.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/lattices.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/matrix_views.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::base::_tags;
using namespace principia::base::_thread_pool;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_lattices;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::quantities::_elementary_functions;

constexpr std::int64_t ε_computation_points = 16;

template<int rows, int columns>
FixedMatrix<cpp_int, rows, columns> ToInt(
    FixedMatrix<cpp_rational, rows, columns> const& m) {
  FixedMatrix<cpp_int, rows, columns> result(uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      auto const& mᵢⱼ = m(i, j);
      DCHECK_EQ(1, denominator(mᵢⱼ));
      result(i, j) = numerator(m(i, j));
    }
  }
  return result;
}

template<int rows, int columns>
FixedMatrix<cpp_rational, rows, columns> ToRational(
    FixedMatrix<cpp_int, rows, columns> const& m) {
  FixedMatrix<cpp_rational, rows, columns> result(uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = m(i, j);
    }
  }
  return result;
}

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

template<std::int64_t zeroes>
cpp_rational GalExhaustiveSearch(std::vector<AccurateFunction> const& functions,
                                 cpp_rational const& starting_argument) {
  CHECK_LT(0, starting_argument);

  // We will look for candidates both above and below |starting_argument|.
  // Note that if |starting_argument| is a power of 2, the increments above
  // and below |starting_argument| are not the same.
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
    if (std::all_of(functions.begin(), functions.end(),
                    [&high_x](AccurateFunction const& f) {
                      return HasDesiredZeroes<zeroes>(f(high_x));
                    })) {
      return high_x;
    }
    high_x += high_increment;
    if (std::all_of(functions.begin(), functions.end(),
                    [&low_x](AccurateFunction const& f) {
                      return HasDesiredZeroes<zeroes>(f(low_x));
                    })) {
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
    cpp_rational const& near_argument,
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
    F[i] = [&functions, i, N, &near_argument](cpp_rational const& t) {
      // Here |t| <= |T|.
      return N * functions[i](near_argument + t / N);
    };
  }
  AccuratePolynomial<cpp_rational, 1> const shift_and_rescale(
      {near_argument, cpp_rational(1, N)});
  for (std::int64_t i = 0; i < polynomials.size(); ++i) {
    P[i] = N * Compose(polynomials[i], shift_and_rescale);
  }

  // Step 2: compute ε.  We don't care too much about its accuracy because in
  // general the largest error is on the boundary of the domain, and anyway ε
  // has virtually no incidence on the value of Mʹ.
  cpp_rational const T_increment = cpp_rational(T, ε_computation_points);
  cpp_bin_float_50 ε = 0;
  for (std::int64_t i = 0; i < 2; ++i) {
    for (cpp_rational t = -T; t <= T; t += T_increment) {
      ε = std::max(ε, abs(F[i](t) - static_cast<cpp_bin_float_50>((*P[i])(t))));
    }
  }
  VLOG(1) << "ε = " << ε;

  // Step 3, first part: compute Mʹ and C.  Give up is C is 0, which may happen
  // if ε is too large.
  auto const Mʹ = static_cast<std::int64_t>(floor(M / (2 + 2 * M * ε)));
  auto const C = 3 * Mʹ;
  if (C == 0) {
    return absl::FailedPreconditionError("Error too large");
  }
  VLOG(1) << "C = " << C;

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
    VLOG(1) << "P̃[" << i << "] = " << *P̃[i];
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
  VLOG(1) << "L = " << L;

  // Step 7: reduce the lattice.
  // The lattice really has integer coefficients, but this is inconvenient to
  // propagate through the matrix algorithms.  (It would require copies instead
  // of views for all the types, not just the ones we use here.)
  Lattice const V = ToInt(LenstraLenstraLovász(ToRational(L)));
  VLOG(1) << "V = " << V;

  // Step 8: find the three shortest vectors of the reduced lattice.  We sort
  // the columns according to the L₂ norm.
  std::array<std::unique_ptr<ColumnView<Lattice const>>, V.columns()> v;
  for (std::int64_t i = 0; i < v.size(); ++i) {
    // TODO(phl): Switch the matrices to use std::int64_t instead of int, and
    // the cast below will go away.
    v[i] = std::make_unique<ColumnView<Lattice const>>(
        ColumnView<Lattice const>{.matrix = V,
                                  .first_row = 0,
                                  .last_row = V.rows() - 1,
                                  .column = static_cast<int>(i)});
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
  FixedMatrix<cpp_rational, 5, dimension> w;
  for (std::int64_t i = 0; i < dimension; ++i) {
    auto const& vᵢ = *v[i];
    VLOG(1) << "v[" << i << "] = " << vᵢ;
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

  if (Q_coefficients[1] == 0) {
      return absl::NotFoundError("No integer zeroes");
  }

  AccuratePolynomial<cpp_rational, 1> const Q({Q_coefficients[0],
                                               Q_coefficients[1]});
  VLOG(1) << "Q = " << Q;

  // Step 11: compute q and find its integer root (singular), if any.
  AccuratePolynomial<cpp_rational, 1> const q =
      Compose(Q, AccuratePolynomial<cpp_rational, 1>({0, cpp_rational(1, T)}));

  cpp_rational const t₀ =
      -std::get<0>(q.coefficients()) / std::get<1>(q.coefficients());
  VLOG(1) << "t₀ = " << t₀;
  if (abs(t₀) > T) {
    return absl::NotFoundError("Out of bounds");
  } else if (denominator(t₀) != 1) {
    return absl::NotFoundError("Noninteger root");
  }

  for (auto const& Fᵢ : F) {
    auto const Fᵢ_t₀ = Fᵢ(t₀);
    auto const Fᵢ_t₀_cmod_1 = Fᵢ_t₀ - round(Fᵢ_t₀);
    if (M * abs(Fᵢ_t₀_cmod_1) >= 1) {
      VLOG(1) << "Fi(t₀) cmod 1 = " << Fᵢ_t₀_cmod_1;
      return absl::NotFoundError("Not enough zeroes");
    }
  }

  return t₀ / N + near_argument;
}

}  // namespace internal
}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
