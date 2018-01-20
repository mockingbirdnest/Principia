#pragma once

#include "numerics/newhall.hpp"

#include <vector>

#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_newhall {

using quantities::Exponentiation;
using quantities::Frequency;
using quantities::Time;

template<typename Vector, typename HeterogeneousCoefficients, int n, int k>
struct Heterogeneizer {
  static void Convert(FixedVector<Vector, n> const& homogeneous_coefficients,
                      Frequency const& scale,
                      Exponentiation<Frequency, k> const& scale_k,
                      HeterogeneousCoefficients& heterogenous_coefficients);
};

template<typename Vector, typename HeterogeneousCoefficients, int n>
struct Heterogeneizer<Vector, HeterogeneousCoefficients, n, n> {
  static void Convert(FixedVector<Vector, n> const& homogeneous_coefficients,
                      Frequency const& time_scale,
                      Exponentiation<Frequency, n> const& scale_n,
                      HeterogeneousCoefficients& heterogenous_coefficients);
};

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Time, degree, Evaluator> Heterogeneize(
    FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale) {
  using P = PolynomialInMonomialBasis<Vector, Time, degree, Evaluator>;
  typename P::Coefficients heterogeneous_coefficients;
  Heterogeneizer<Vector, typename P::Coefficients, degree + 1, 0>::Convert(
      homogeneous_coefficients,
      scale,
      /*scale_0=*/1.0,
      heterogeneous_coefficients);
  return P(heterogeneous_coefficients);
}

template<typename Vector, typename HeterogeneousCoefficients, int n, int k>
void Heterogeneizer<Vector, HeterogeneousCoefficients, n, k>::Convert(
    FixedVector<Vector, n> const& homogeneous_coefficients,
    Frequency const& scale,
    Exponentiation<Frequency, k> const& scale_k,
    HeterogeneousCoefficients& heterogenous_coefficients) {
  std::get<k>(heterogenous_coefficients) =
      homogeneous_coefficients[k] * scale_k;
  Heterogeneizer<Vector, HeterogeneousCoefficients, n, k + 1>::Convert(
      homogeneous_coefficients,
      scale,
      scale_k * scale,
      heterogenous_coefficients);
}

template<typename Vector, typename HeterogeneousCoefficients, int n>
void Heterogeneizer<Vector, HeterogeneousCoefficients, n, n>::Convert(
    FixedVector<Vector, n> const& homogeneous_coefficients,
    Frequency const& scale,
    Exponentiation<Frequency, n> const& scale_n,
    HeterogeneousCoefficients& heterogenous_coefficients) {
  std::get<n>(heterogenous_coefficients) =
      homogeneous_coefficients[n] * scale_n;
}

template<typename Vector>
ЧебышёвSeries<Vector>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<Vector> const& q,
                                   std::vector<Variation<Vector>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Vector& error_estimate) {
  // Only supports 8 divisions for now.
  int const divisions = 8;
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Vector, 2 * divisions + 2> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i];
    qv[j + 1] = v[i] * duration_over_two;
  }

  std::vector<Vector> coefficients;
  coefficients.reserve(degree);
  switch (degree) {
    case 3:
      coefficients = newhall_c_matrix_чебышёв_degree_3_divisions_8_w04 * qv;
      break;
    case 4:
      coefficients = newhall_c_matrix_чебышёв_degree_4_divisions_8_w04 * qv;
      break;
    case 5:
      coefficients = newhall_c_matrix_чебышёв_degree_5_divisions_8_w04 * qv;
      break;
    case 6:
      coefficients = newhall_c_matrix_чебышёв_degree_6_divisions_8_w04 * qv;
      break;
    case 7:
      coefficients = newhall_c_matrix_чебышёв_degree_7_divisions_8_w04 * qv;
      break;
    case 8:
      coefficients = newhall_c_matrix_чебышёв_degree_8_divisions_8_w04 * qv;
      break;
    case 9:
      coefficients = newhall_c_matrix_чебышёв_degree_9_divisions_8_w04 * qv;
      break;
    case 10:
      coefficients = newhall_c_matrix_чебышёв_degree_10_divisions_8_w04 * qv;
      break;
    case 11:
      coefficients = newhall_c_matrix_чебышёв_degree_11_divisions_8_w04 * qv;
      break;
    case 12:
      coefficients = newhall_c_matrix_чебышёв_degree_12_divisions_8_w04 * qv;
      break;
    case 13:
      coefficients = newhall_c_matrix_чебышёв_degree_13_divisions_8_w04 * qv;
      break;
    case 14:
      coefficients = newhall_c_matrix_чебышёв_degree_14_divisions_8_w04 * qv;
      break;
    case 15:
      coefficients = newhall_c_matrix_чебышёв_degree_15_divisions_8_w04 * qv;
      break;
    case 16:
      coefficients = newhall_c_matrix_чебышёв_degree_16_divisions_8_w04 * qv;
      break;
    case 17:
      coefficients = newhall_c_matrix_чебышёв_degree_17_divisions_8_w04 * qv;
      break;
    default:
      LOG(FATAL) << "Unexpected degree " << degree;
      break;
  }
  CHECK_EQ(degree + 1, coefficients.size());
  error_estimate = coefficients[degree];
  return ЧебышёвSeries<Vector>(coefficients, t_min, t_max);
}

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Time, degree, Evaluator>
NewhallApproximationInMonomialBasis(std::vector<Vector> const& q,
                                    std::vector<Variation<Vector>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Vector& error_estimate) {
  // Only supports 8 divisions for now.
  int const divisions = 8;
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Vector, 2 * divisions + 2> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i];
    qv[j + 1] = v[i] * duration_over_two;
  }
  Frequency const one_over_duration_over_two = 1.0 / duration_over_two;

  switch (degree) {
    case 3:
      return Heterogeneize<Vector, 3, Evaluator>(
          newhall_c_matrix_monomial_degree_3_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 4:
      return Heterogeneize<Vector, 4, Evaluator>(
          newhall_c_matrix_monomial_degree_4_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 5:
      return Heterogeneize<Vector, 5, Evaluator>(
          newhall_c_matrix_monomial_degree_5_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 6:
      return Heterogeneize<Vector, 6, Evaluator>(
          newhall_c_matrix_monomial_degree_6_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 7:
      return Heterogeneize<Vector, 7, Evaluator>(
          newhall_c_matrix_monomial_degree_7_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 8:
      return Heterogeneize<Vector, 8, Evaluator>(
          newhall_c_matrix_monomial_degree_8_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 9:
      return Heterogeneize<Vector, 9, Evaluator>(
          newhall_c_matrix_monomial_degree_9_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 10:
      return Heterogeneize<Vector, 10, Evaluator>(
          newhall_c_matrix_monomial_degree_10_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 11:
      return Heterogeneize<Vector, 11, Evaluator>(
          newhall_c_matrix_monomial_degree_11_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 12:
      return Heterogeneize<Vector, 12, Evaluator>(
          newhall_c_matrix_monomial_degree_12_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 13:
      return Heterogeneize<Vector, 13, Evaluator>(
          newhall_c_matrix_monomial_degree_13_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 14:
      return Heterogeneize<Vector, 14, Evaluator>(
          newhall_c_matrix_monomial_degree_14_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 15:
      return Heterogeneize<Vector, 15, Evaluator>(
          newhall_c_matrix_monomial_degree_15_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 16:
      return Heterogeneize<Vector, 16, Evaluator>(
          newhall_c_matrix_monomial_degree_16_divisions_8_w04 * qv,
          one_over_duration_over_two);
    case 17:
      return Heterogeneize<Vector, 17, Evaluator>(
          newhall_c_matrix_monomial_degree_17_divisions_8_w04 * qv,
          one_over_duration_over_two);
    default:
      LOG(FATAL) << "Unexpected degree " << degree;
      break;
  }
  //error estimate
}

}  // namespace internal_newhall
}  // namespace numerics
}  // namespace principia
