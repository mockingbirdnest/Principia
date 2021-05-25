#pragma once

#include "numerics/newhall.hpp"

#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_newhall {

using base::make_not_null_unique;
using geometry::Barycentre;
using quantities::Exponentiation;
using quantities::Frequency;
using quantities::Time;

// Only supports 8 divisions for now.
constexpr int divisions = 8;

template<typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Argument, Instant, degree, Evaluator> Dehomogeneize(
    FixedVector<Argument, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin);

template<typename Argument,
         typename DehomogeneizedCoefficients, int degree, int k>
struct Dehomogeneizer {
  static void Convert(
      FixedVector<Argument, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, k> const& scale_k,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Argument,
         typename DehomogeneizedCoefficients, int degree>
struct Dehomogeneizer<Argument, DehomogeneizedCoefficients, degree, degree> {
  static void Convert(
      FixedVector<Argument, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, degree> const& scale_degree,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Argument, Instant, degree, Evaluator> Dehomogeneize(
    FixedVector<Argument, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin) {
  using P = PolynomialInMonomialBasis<Argument, Instant, degree, Evaluator>;
  typename P::Coefficients dehomogeneized_coefficients;
  Dehomogeneizer<Argument, typename P::Coefficients, degree, /*k=*/0>::Convert(
      homogeneous_coefficients,
      scale,
      /*scale_k=*/1.0,
      dehomogeneized_coefficients);
  return P(dehomogeneized_coefficients, origin);
}

template<typename Argument,
         typename DehomogeneizedCoefficients, int degree, int k>
void Dehomogeneizer<Argument, DehomogeneizedCoefficients, degree, k>::Convert(
    FixedVector<Argument, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Exponentiation<Frequency, k> const& scale_k,
    DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<k>(dehomogeneized_coefficients) =
      homogeneous_coefficients[k] * scale_k;
  Dehomogeneizer<Argument, DehomogeneizedCoefficients, degree, k + 1>::Convert(
      homogeneous_coefficients,
      scale,
      scale_k * scale,
      dehomogeneized_coefficients);
}

template<typename Argument,
         typename DehomogeneizedCoefficients, int degree>
void Dehomogeneizer<Argument, DehomogeneizedCoefficients, degree, degree>::
Convert(FixedVector<Argument, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        Exponentiation<Frequency, degree> const& scale_degree,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<degree>(dehomogeneized_coefficients) =
      homogeneous_coefficients[degree] * scale_degree;
}

template<typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator {
  static FixedVector<Argument, degree + 1> HomogeneousCoefficients(
      FixedVector<Argument, 2 * divisions + 2> const& qv,
      Argument& error_estimate);
};

#define PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(degree)                  \
  template<typename Argument,                                                  \
           template<typename, typename, int>                                   \
           class Evaluator>                                                    \
  struct NewhallAppromixator<Argument, (degree), Evaluator> {                  \
    static FixedVector<Argument, ((degree) + 1)> HomogeneousCoefficients(      \
        FixedVector<Argument, 2 * divisions + 2> const& qv,                    \
        Argument& error_estimate) {                                            \
      error_estimate =                                                         \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04           \
              .row<(degree)>() *                                               \
          qv;                                                                  \
      return newhall_c_matrix_monomial_degree_##degree##_divisions_8_w04 * qv; \
    }                                                                          \
  }

PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(3);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(4);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(5);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(6);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(7);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(8);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(9);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(10);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(11);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(12);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(13);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(14);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(15);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(16);
PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(17);

#undef PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(degree)     \
  case (degree):                                                          \
    coefficients = std::vector<Argument>(                                 \
        newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04 * qv); \
    break

template<typename Argument>
ЧебышёвSeries<Argument>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<Argument> const& q,
                                   std::vector<Variation<Argument>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Difference<Argument>& error_estimate) {
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Argument const origin{};
  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Argument, 2 * divisions + 2> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i] - origin;
    qv[j + 1] = v[i] * duration_over_two;
  }

  std::vector<Argument> coefficients;
  coefficients.reserve(degree);
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
  CHECK_EQ(degree + 1, coefficients.size());
  error_estimate = coefficients[degree];
  return ЧебышёвSeries<Difference<Argument>>(coefficients, t_min, t_max);
}

#undef PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE

template<typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Argument, Instant, degree, Evaluator>
NewhallApproximationInMonomialBasis(std::vector<Argument> const& q,
                                    std::vector<Variation<Argument>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Difference<Argument>& error_estimate) {
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Argument const origin{};
  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Difference<Argument>, 2 * divisions + 2> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i] - origin;
    qv[j + 1] = v[i] * duration_over_two;
  }

  Instant const t_mid = Barycentre<Instant, double>({t_min, t_max}, {1, 1});
  return Dehomogeneize<Difference<Argument>, degree, Evaluator>(
             NewhallAppromixator<Difference<Argument>, degree, Evaluator>::
                 HomogeneousCoefficients(qv, error_estimate),
             /*scale=*/1.0 / duration_over_two,
             t_mid);
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(degree)      \
  case (degree):                                                            \
    return make_not_null_unique<                                            \
        PolynomialInMonomialBasis<Argument, Instant, (degree), Evaluator>>( \
        NewhallApproximationInMonomialBasis<Argument, (degree), Evaluator>( \
            q, v, t_min, t_max, error_estimate))

template<typename Argument, template<typename, typename, int> class Evaluator>
not_null<std::unique_ptr<Polynomial<Argument, Instant>>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<Argument> const& q,
                                    std::vector<Variation<Argument>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Difference<Argument>& error_estimate) {
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

}  // namespace internal_newhall
}  // namespace numerics
}  // namespace principia
