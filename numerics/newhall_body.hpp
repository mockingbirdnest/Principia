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

using geometry::Barycentre;
using quantities::Exponentiation;
using quantities::Frequency;
using quantities::Time;
using namespace principia::base::_not_null;

// Only supports 8 divisions for now.
constexpr int divisions = 8;

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator> Dehomogeneize(
    FixedVector<Value, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin);

template<typename Value,
         typename DehomogeneizedCoefficients, int degree, int k>
struct Dehomogeneizer {
  static void Convert(
      FixedVector<Value, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, k> const& scale_k,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Value,
         typename DehomogeneizedCoefficients, int degree>
struct Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, degree> {
  static void Convert(
      FixedVector<Value, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      Exponentiation<Frequency, degree> const& scale_degree,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator> Dehomogeneize(
    FixedVector<Value, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    Instant const& origin) {
  using P = PolynomialInMonomialBasis<Value, Instant, degree, Evaluator>;
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
void Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, k>::Convert(
    FixedVector<Value, degree + 1> const& homogeneous_coefficients,
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
Convert(FixedVector<Value, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        Exponentiation<Frequency, degree> const& scale_degree,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<degree>(dehomogeneized_coefficients) =
      homogeneous_coefficients[degree] * scale_degree;
}

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator {
  static FixedVector<Value, degree + 1> HomogeneousCoefficients(
      FixedVector<Value, 2 * divisions + 2> const& qv,
      Value& error_estimate);
};

#define PRINCIPIA_NEWHALL_APPROXIMATOR_SPECIALIZATION(degree)                  \
  template<typename Value, template<typename, typename, int> class Evaluator>  \
  struct NewhallAppromixator<Value, (degree), Evaluator> {                     \
    static FixedVector<Value, ((degree) + 1)> HomogeneousCoefficients(         \
        FixedVector<Value, 2 * divisions + 2> const& qv,                       \
        Value& error_estimate) {                                               \
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
    coefficients = std::vector<Vector>(                                   \
        newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04 * qv); \
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

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator>
NewhallApproximationInMonomialBasis(std::vector<Value> const& q,
                                    std::vector<Variation<Value>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Difference<Value>& error_estimate) {
  CHECK_EQ(divisions + 1, q.size());
  CHECK_EQ(divisions + 1, v.size());

  Value const origin{};
  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Difference<Value>, 2 * divisions + 2> qv;
  for (int i = 0, j = 2 * divisions;
       i < divisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i] - origin;
    qv[j + 1] = v[i] * duration_over_two;
  }

  Instant const t_mid = Barycentre<Instant, double>({t_min, t_max}, {1, 1});
  return origin +
         Dehomogeneize<Difference<Value>, degree, Evaluator>(
             NewhallAppromixator<Difference<Value>, degree, Evaluator>::
                 HomogeneousCoefficients(qv, error_estimate),
             /*scale=*/1.0 / duration_over_two,
             t_mid);
}

#define PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE(degree)   \
  case (degree):                                                         \
    return make_not_null_unique<                                         \
        PolynomialInMonomialBasis<Value, Instant, (degree), Evaluator>>( \
        NewhallApproximationInMonomialBasis<Value, (degree), Evaluator>( \
            q, v, t_min, t_max, error_estimate))

template<typename Value, template<typename, typename, int> class Evaluator>
not_null<std::unique_ptr<Polynomial<Value, Instant>>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<Value> const& q,
                                    std::vector<Variation<Value>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
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

}  // namespace internal_newhall
}  // namespace numerics
}  // namespace principia
