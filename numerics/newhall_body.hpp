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

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Instant, degree, Evaluator> Dehomogeneize(
    FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin);

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree, int k>
struct Dehomogeneizer {
  static void Convert(
      FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, k> const& scale_k,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree>
struct Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, degree> {
  static void Convert(
      FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, degree> const& scale_degree,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Instant, degree, Evaluator> Dehomogeneize(
    FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin) {
  using P = PolynomialInMonomialBasis<Vector, Instant, degree, Evaluator>;
  typename P::Coefficients dehomogeneized_coefficients;
  Dehomogeneizer<Vector, typename P::Coefficients, degree, /*k=*/0>::Convert(
      homogeneous_coefficients,
      scale,
      /*scale_k=*/1.0,
      dehomogeneized_coefficients);
  return P(dehomogeneized_coefficients, origin);
}

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree, int k>
void Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, k>::Convert(
    FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Exponentiation<Frequency, k> const& scale_k,
    DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<k>(dehomogeneized_coefficients) =
      homogeneous_coefficients[k] * scale_k;
  Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, k + 1>::Convert(
      homogeneous_coefficients,
      scale,
      scale_k * scale,
      dehomogeneized_coefficients);
}

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree>
void Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, degree>::
Convert(FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        Exponentiation<Frequency, degree> const& scale_degree,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<degree>(dehomogeneized_coefficients) =
      homogeneous_coefficients[degree] * scale_degree;
}

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator {
  static FixedVector<Vector, degree + 1> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv,
      Vector& error_estimate);
};

#define PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(degree)                  \
  template<typename Vector,                                                    \
           template<typename, typename, int> class Evaluator>                  \
  struct NewhallAppromixator<Vector, (degree), Evaluator> {                    \
    static FixedVector<Vector, ((degree) + 1)> HomogeneousCoefficients(        \
        FixedVector<Vector, 2 * divisions + 2> const& qv,                      \
        Vector& error_estimate) {                                              \
      error_estimate =                                                         \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04.row<      \
              (degree)>() * qv;                                                \
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

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE(degree)    \
  case (degree):                                                         \
    coefficients =                                                       \
        newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04 * qv; \
    break

template<typename Vector>
ЧебышёвSeries<Vector>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<Vector> const& q,
                                   std::vector<Variation<Vector>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Vector& error_estimate) {
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
  return ЧебышёвSeries<Vector>(coefficients, t_min, t_max);
}

#undef PRINCIPIA_NEWHALL_APPROXIMATION_IN_ЧЕБЫШЁВ_BASIS_CASE

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Instant, degree, Evaluator>
NewhallApproximationInMonomialBasis(std::vector<Vector> const& q,
                                    std::vector<Variation<Vector>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Vector& error_estimate) {
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

  Instant const t_mid = Barycentre<Instant, double>({t_min, t_max}, {1, 1});
  return Dehomogeneize<Vector, degree, Evaluator>(
             NewhallAppromixator<Vector, degree, Evaluator>::
                 HomogeneousCoefficients(qv, error_estimate),
             /*scale=*/1.0 / duration_over_two,
             t_mid);
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(degree)    \
  case (degree):                                                          \
    return make_not_null_unique<                                          \
        PolynomialInMonomialBasis<Vector, Instant, (degree), Evaluator>>( \
        NewhallApproximationInMonomialBasis<Vector, (degree), Evaluator>( \
            q, v,                                                         \
            t_min, t_max,                                                 \
            error_estimate))

template<typename Vector, template<typename, typename, int> class Evaluator>
not_null<std::unique_ptr<Polynomial<Vector, Instant>>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<Vector> const& q,
                                    std::vector<Variation<Vector>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Vector& error_estimate) {
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
