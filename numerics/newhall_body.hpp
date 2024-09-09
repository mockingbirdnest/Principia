#pragma once

#include "numerics/newhall.hpp"

#include <memory>
#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/newhall_matrices.mathematica.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _newhall {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_newhall_matrices;
using namespace principia::quantities::_quantities;

// Only supports 8 divisions for now.
constexpr int divisions = 8;

// Logically `QV` should be heterogeneous, as it contains positions (or
// displacements) and velocities.  However, this would require giving dimensions
// to the derivatives of the Чебышёв polynomials.  Let's no go there, let's do a
// bit of type decay instead.
template<typename Value>
using QV = std::array<Value, 2 * divisions + 2>;

// A helper to unroll the dot product between an array-like object (which must
// have an operator[]) and a `QV`.
template<int index = 2 * divisions + 1>
struct DotProduct {
  template<typename Left, typename RightElement>
  static RightElement Compute(Left const& left, QV<RightElement> const& right);
};

template<>
struct DotProduct<0> {
  template<typename Left, typename RightElement>
  static RightElement Compute(Left const& left, QV<RightElement> const& right);
};

template<int index>
template<typename Left, typename RightElement>
RightElement DotProduct<index>::Compute(Left const& left,
                                        QV<RightElement> const& right) {
  return left[index] * right[index] +
         DotProduct<index - 1>::Compute(left, right);
}

template<typename Left, typename RightElement>
RightElement DotProduct<0>::Compute(Left const& left,
                                    QV<RightElement> const& right) {
  return left[0] * right[0];
}

// Fills `result` (which must be array-like) with the result of multiplying a
// matrix with a `QV`.
template<int degree, typename RightElement, typename Result>
void Multiply(FixedMatrix<double, degree + 1, 2 * divisions + 2> const& left,
              QV<RightElement> const& right,
              Result& result) {
  auto const* row = left.template row<0>();
  for (int i = 0; i < degree + 1; ++i) {
    result[i] = DotProduct<>::Compute(row, right);
    row += 2 * divisions + 2;
  }
}

template<typename Value, int degree>
PolynomialInMonomialBasis<Value, Instant, degree> Dehomogeneize(
    std::array<Value, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin);

template<typename Value,
         typename DehomogeneizedCoefficients, int degree, int k>
struct Dehomogeneizer {
  static void Convert(
      std::array<Value, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, k> const& scale_k,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Value,
         typename DehomogeneizedCoefficients, int degree>
struct Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, degree> {
  static void Convert(
      std::array<Value, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, degree> const& scale_degree,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Value, int degree>
PolynomialInMonomialBasis<Value, Instant, degree> Dehomogeneize(
    std::array<Value, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin) {
  using P = PolynomialInMonomialBasis<Value, Instant, degree>;
  typename P::Coefficients dehomogeneized_coefficients;
  Dehomogeneizer<Value, typename P::Coefficients, degree, /*k=*/0>::Convert(
      homogeneous_coefficients,
      scale,
      /*scale_k=*/1.0,
      dehomogeneized_coefficients);
  return P(dehomogeneized_coefficients, origin);
}

template<typename Value,
         typename DehomogeneizedCoefficients, int degree, int k>
void Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, k>::
Convert(std::array<Value, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        Exponentiation<Frequency, k> const& scale_k,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<k>(dehomogeneized_coefficients) =
      homogeneous_coefficients[k] * scale_k;
  Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, k + 1>::Convert(
      homogeneous_coefficients,
      scale,
      scale_k * scale,
      dehomogeneized_coefficients);
}

template<typename Value,
         typename DehomogeneizedCoefficients, int degree>
void Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, degree>::
Convert(std::array<Value, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        Exponentiation<Frequency, degree> const& scale_degree,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<degree>(dehomogeneized_coefficients) =
      homogeneous_coefficients[degree] * scale_degree;
}

template<typename Value, int degree>
struct NewhallЧебышёвApproximator {
  static std::array<Value, degree + 1> HomogeneousCoefficients(
      QV<Value> const& qv,
      Value& error_estimate);
};

template<typename Value, int degree>
struct NewhallMonomialApproximator {
  static std::array<Value, degree + 1> HomogeneousCoefficients(
      QV<Value> const& qv,
      Value& error_estimate);
};

#define PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(degree) \
  template<typename Value>                                            \
  struct NewhallЧебышёвApproximator<Value, (degree)> {                \
    static std::array<Value, (degree) + 1> HomogeneousCoefficients(   \
        QV<Value> const& qv,                                          \
        Value& error_estimate) {                                      \
      std::array<Value, degree + 1> result;                           \
      Multiply<(degree)>(                                             \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04, \
          qv,                                                         \
          result);                                                    \
      error_estimate = result[degree];                                \
      return result;                                                  \
    }                                                                 \
  }

// The error estimate must be computed in the Чебышёв basis because the elements
// of the monomial basis are not bounded.
#define PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(degree) \
  template<typename Value>                                             \
  struct NewhallMonomialApproximator<Value, (degree)> {                \
    static std::array<Value, (degree) + 1> HomogeneousCoefficients(    \
        QV<Value> const& qv,                                           \
        Value& error_estimate) {                                       \
      error_estimate = DotProduct<>::Compute(                          \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04   \
              .row<(degree)>(),                                        \
          qv);                                                         \
      std::array<Value, (degree) + 1> result;                          \
      Multiply<(degree)>(                                              \
          newhall_c_matrix_monomial_degree_##degree##_divisions_8_w04, \
          qv,                                                          \
          result);                                                     \
      return result;                                                   \
    }                                                                  \
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
NewhallApproximationInЧебышёвBasis(std::vector<Value> const& q,
                                   std::vector<Variation<Value>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Value& error_estimate) {
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  QV<Value> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i];
    qv[j + 1] = v[i] * duration_over_two;
  }

  auto const coefficients =
      NewhallЧебышёвApproximator<Difference<Value>, degree>::
          HomogeneousCoefficients(qv, error_estimate);
  return PolynomialInЧебышёвBasis<Value, Instant, degree>(
      coefficients, t_min, t_max);
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(degree) \
  case (degree):                                                      \
    return make_not_null_unique<                                      \
        PolynomialInЧебышёвBasis<Value, Instant, (degree)>>(          \
        NewhallApproximationInЧебышёвBasis<Value, (degree)>(          \
            q, v, t_min, t_max, error_estimate))

template<typename Value>
not_null<std::unique_ptr<PolynomialInЧебышёвBasis<Value, Instant>>>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<Value> const& q,
                                   std::vector<Variation<Value>> const& v,
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

template<typename Value, int degree>
PolynomialInMonomialBasis<Value, Instant, degree>
NewhallApproximationInMonomialBasis(std::vector<Value> const& q,
                                    std::vector<Variation<Value>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Policy const& policy,
                                    Difference<Value>& error_estimate) {
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Value const origin{};
  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  QV<Difference<Value>> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i] - origin;
    qv[j + 1] = v[i] * duration_over_two;
  }

  Instant const t_mid = Barycentre({t_min, t_max});
  return policy.WithEvaluator(
      origin +
      Dehomogeneize<Difference<Value>, degree>(
          NewhallMonomialApproximator<Difference<Value>, degree>::
              HomogeneousCoefficients(qv, error_estimate),
          /*scale=*/1.0 / duration_over_two,
          t_mid));
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(degree) \
  case (degree):                                                       \
    return make_not_null_unique<                                       \
        PolynomialInMonomialBasis<Value, Instant, (degree)>>(          \
        NewhallApproximationInMonomialBasis<Value, (degree)>(          \
            q, v, t_min, t_max, policy, error_estimate))

template<typename Value>
not_null<std::unique_ptr<Polynomial<Value, Instant>>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<Value> const& q,
                                    std::vector<Variation<Value>> const& v,
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
      break;
  }
}

#undef PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE

}  // namespace internal
}  // namespace _newhall
}  // namespace numerics
}  // namespace principia
