#pragma once

#include "numerics/newhall.hpp"

#include <memory>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/macros.hpp"  // 🧙 For FORCE_INLINE.
#include "base/tags.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "numerics/elementary_functions.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/newhall_matrices.mathematica.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _newhall {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::base::_tags;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_newhall_matrices;
using namespace principia::quantities::_quantities;

// Only supports 8 divisions for now.
constexpr int divisions = 8;

// The input is given as typed pairs of type `QV` with an argument in
// [t_min, t_max].  However, the argument must be rescaled to [-1, 1] because
// that's the interval over which the Чебышёв polynomials are defined ([New89],
// p. 305).  After rescaling, the typed pairs contain two elements of type
// `Value`.
template<typename Value>
using RescaledQV = DirectSum<Value, Value>;

template<typename Value>
using RescaledQVs = std::array<RescaledQV<Value>, divisions + 1>;

using NewhallMatrixElement = DirectSum<double, double>;

// Multiplies a row of a matrix (given as a pointer to its first element) by a
// column vector.  In [New89], p. 308, the former is a row of `C₁⁻¹C₂` , the
// latter is the vector `f`.
template<typename Value>
FORCE_INLINE constexpr Value MultiplyMatrixRowByColumnVector(
    NewhallMatrixElement const* const left,
    RescaledQVs<Value> const& right) {
  Value result{};

  FORCE_INLINE_CALLS
  for_integer_range<0, divisions + 1>::loop([&]<int i> FORCE_INLINE {
    auto const& [l0, l1] = left[i];
    auto const& [r0, r1] = right[i];
    // This computation preserves the accuracy obtained with the previous
    // implementation.
    result += l0 * r0;
    result += l1 * r1;
  });
  return result;
}

// Similar to the previous function, but performs the multiplication for all
// rows of the matrix, returning the result as an array.
template<int degree, typename Value>
FORCE_INLINE constexpr
std::array<Value, degree + 1> MultiplyMatrixByColumnVector(
    FixedMatrix<NewhallMatrixElement, degree + 1, divisions + 1> const& left,
    RescaledQVs<Value> const& right) {
  std::array<Value, degree + 1> result;
  for (std::int64_t i = 0; i <= degree; ++i) {
    auto const* row = left.row(i);
    result[i] = MultiplyMatrixRowByColumnVector(row, right);
  }
  return result;
}

// Given the coefficients of a homogeneous polynomial with argument `double`,
// returns the coefficients adjusted to the argument interval [t_min, t_max]
// (specified by `scale`) and the origin of the approximation.  It is necessary
// to go through a homegeneous polynomial because the tables cannot account for
// the width of the interval [t_min, t_max].
template<typename Value, int degree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE constexpr
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator> Dehomogeneize(
    std::array<Difference<Value>, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Value const& value_origin,
    Instant const& argument_origin) {
  using P = PolynomialInMonomialBasis<Value, Instant, degree, Evaluator>;
  typename P::Coefficients dehomogeneized_coefficients(uninitialized);

  FORCE_INLINE_CALLS
  for_integer_range<0, degree + 1>::loop([&]<int k> FORCE_INLINE {
    if constexpr (k == 0) {
      get<k>(dehomogeneized_coefficients) =
          homogeneous_coefficients[k] + value_origin;
    } else {
      get<k>(dehomogeneized_coefficients) =
          homogeneous_coefficients[k] * Pow<k>(scale);
    }
  });

  return P(dehomogeneized_coefficients, argument_origin);
}

template<typename Value, int degree>
struct NewhallЧебышёвApproximator {
  static std::array<Value, degree + 1> ComputeCoefficients(
      RescaledQVs<Value> const& rqvs,
      Value& error_estimate);
};

template<typename Value, int degree>
struct NewhallMonomialApproximator {
  static std::array<Value, degree + 1> ComputeHomogeneousCoefficients(
      RescaledQVs<Value> const& rqvs,
      Value& error_estimate);
};

#define PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(degree)        \
  template<typename Value>                                                   \
  struct NewhallЧебышёвApproximator<Value, (degree)> {                       \
    FORCE_INLINE static constexpr std::array<Value, (degree) + 1>            \
    ComputeCoefficients(RescaledQVs<Value> const& rqvs,                      \
                        Value& error_estimate) {                             \
      auto const result = MultiplyMatrixByColumnVector<(degree)>(            \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04, rqvs); \
      error_estimate = result[(degree)];                                     \
      return result;                                                         \
    }                                                                        \
  }

// The error estimate must be computed in the Чебышёв basis because the elements
// of the monomial basis are not bounded.
#define PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(degree)        \
  template<typename Value>                                                    \
  struct NewhallMonomialApproximator<Value, (degree)> {                       \
    FORCE_INLINE static constexpr std::array<Value, (degree) + 1>             \
    ComputeHomogeneousCoefficients(RescaledQVs<Value> const& rqvs,            \
                                   Value& error_estimate) {                   \
      error_estimate = MultiplyMatrixRowByColumnVector(                       \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04.row(     \
              degree),                                                        \
          rqvs);                                                              \
      return MultiplyMatrixByColumnVector<(degree)>(                          \
          newhall_c_matrix_monomial_degree_##degree##_divisions_8_w04, rqvs); \
    }                                                                         \
  }

PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(3);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(4);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(5);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(6);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(7);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(8);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(9);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(10);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(11);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(12);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(13);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(14);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(15);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(16);
PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(17);

PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(3);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(4);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(5);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(6);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(7);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(8);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(9);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(10);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(11);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(12);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(13);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(14);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(15);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(16);
PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(17);

#undef PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION
#undef PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION

template<typename Value, int degree>
PolynomialInЧебышёвBasis<Value, Instant, degree>
NewhallApproximationInЧебышёвBasis(std::vector<QV<Value>> const& qvs,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Value& error_estimate) {
  CHECK_EQ(divisions + 1, qvs.size());

  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  RescaledQVs<Difference<Value>> rqvs;
  for (int i = 0, j = divisions; i < divisions + 1 && j >= 0; ++i, --j) {
    auto const& [q, v] = qvs[i];
    rqvs[j] = {q, v * duration_over_two};
  }

  auto const coefficients =
      NewhallЧебышёвApproximator<
          Difference<Value>, degree>::ComputeCoefficients(rqvs, error_estimate);
  return PolynomialInЧебышёвBasis<Value, Instant, degree>(
      coefficients, t_min, t_max);
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(degree) \
  case (degree):                                                      \
    return make_not_null_unique<                                      \
        PolynomialInЧебышёвBasis<Value, Instant, (degree)>>(          \
        NewhallApproximationInЧебышёвBasis<Value, (degree)>(          \
            qvs, t_min, t_max, error_estimate))

template<typename Value>
not_null<std::unique_ptr<PolynomialInЧебышёвBasis<Value, Instant>>>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<QV<Value>> const& qvs,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Value& error_estimate) {
  switch (degree) {
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(3);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(4);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(5);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(6);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(7);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(8);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(9);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(10);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(11);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(12);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(13);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(14);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(15);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(16);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(17);
    default:
      LOG(FATAL) << "Unexpected degree " << degree;
      break;
  }
}

#undef PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE

template<typename Value, int degree,
         template<typename, typename, int> typename Evaluator>
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator>
NewhallApproximationInMonomialBasis(std::vector<QV<Value>> const& qvs,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Difference<Value>& error_estimate) {
  CHECK_EQ(divisions + 1, qvs.size());

  Value const origin{};
  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  RescaledQVs<Difference<Value>> rqvs;
  for (int i = 0, j = divisions;
       i < divisions + 1 && j >= 0;
       ++i, --j) {
    auto const& [q, v] = qvs[i];
    rqvs[j] = {q - origin, v * duration_over_two};
  }

  Instant const t_mid = Barycentre({t_min, t_max});
  return Dehomogeneize<Value, degree, Evaluator>(
      NewhallMonomialApproximator<Difference<Value>, degree>::
          ComputeHomogeneousCoefficients(rqvs, error_estimate),
      /*scale=*/1.0 / duration_over_two,
      /*value_origin=*/origin,
      /*argument_origin=*/t_mid);
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(degree) \
  case (degree):                                                       \
    return policy.WithEvaluator(                                       \
        NewhallApproximationInMonomialBasis<Value,                     \
                                            (degree),                  \
                                            DefaultEvaluator>(         \
            qvs, t_min, t_max, error_estimate))

template<typename Value>
not_null<std::unique_ptr<Polynomial<Value, Instant>>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<QV<Value>> const& qvs,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Policy const& policy,
                                    Difference<Value>& error_estimate) {
  switch (degree) {
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(3);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(4);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(5);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(6);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(7);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(8);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(9);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(10);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(11);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(12);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(13);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(14);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(15);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(16);
    PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(17);
    default:
      LOG(FATAL) << "Unexpected degree " << degree;
#if PRINCIPIA_COMPILER_MSVC && (_MSC_FULL_VER == 194'435'224)
      std::abort();
#endif
  }
}

#undef PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE

}  // namespace internal
}  // namespace _newhall
}  // namespace numerics
}  // namespace principia
