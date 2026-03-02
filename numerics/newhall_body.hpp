#pragma once

#include "numerics/newhall.hpp"

#include <memory>
#include <vector>

#include "base/for_all_of.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/newhall_matrices.mathematica.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _newhall {
namespace internal {

using namespace principia::base::_for_all_of;
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
//TODO(phl)comment
template<typename Value>
using HomogeneousQV = DirectSum<Value, Value>;

template<typename Value>
using HomogeneousQVs = std::array<HomogeneousQV<Value>, divisions + 1>;

// TODO(phl)comment,name
template<typename Value>
Value DotProduct(DirectSum<double, double> const* const left,
                 HomogeneousQVs<Value> const& right) {
  Value result{};
  for_integer_range<0, divisions + 1>(
      [&]<int i> { result += left[i] * right[i]; });
  return result;
}

// TODO(phl)comment,name
template<int degree, typename Value>
std::array<Value, degree + 1> Multiply(
    FixedMatrix<DirectSum<double, double>, degree + 1, divisions + 1> const&
        left,
    HomogeneousQVs<Value> const& right) {
  std::array<Value, degree + 1> result;
  for_integer_range<0, degree + 1>([&]<int i> {
    // TODO(phl)not a ptr
    auto const* row = left.template row<i>();
    result[i] = DotProduct(row, right);
  });
  return result;
}

template<typename Value, int degree,
         template<typename, typename, int> typename Evaluator>
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator> Dehomogeneize(
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

template<typename Value, int degree,
         template<typename, typename, int> typename Evaluator>
PolynomialInMonomialBasis<Value, Instant, degree, Evaluator> Dehomogeneize(
    std::array<Value, degree + 1> const& homogeneous_coefficients,
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
void Dehomogeneizer<Value, DehomogeneizedCoefficients, degree, k>::
Convert(std::array<Value, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        Exponentiation<Frequency, k> const& scale_k,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  get<k>(dehomogeneized_coefficients) =
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
  get<degree>(dehomogeneized_coefficients) =
      homogeneous_coefficients[degree] * scale_degree;
}

template<typename Value, int degree>
struct NewhallЧебышёвApproximator {
  static std::array<Value, degree + 1> HomogeneousCoefficients(
      HomogeneousQVs<Value> const& hqvs,
      Value& error_estimate);
};

template<typename Value, int degree>
struct NewhallMonomialApproximator {
  static std::array<Value, degree + 1> HomogeneousCoefficients(
      HomogeneousQVs<Value> const& qvs,
      Value& error_estimate);
};

#define PRINCIPIA_NEWHALL_ЧЕБЫШЁВ_APPROXIMATOR_SPECIALIZATION(degree)        \
  template<typename Value>                                                   \
  struct NewhallЧебышёвApproximator<Value, (degree)> {                       \
    static std::array<Value, (degree) + 1> HomogeneousCoefficients(          \
        HomogeneousQVs<Value> const& hqvs,                                   \
        Value& error_estimate) {                                             \
      auto const result = Multiply<(degree)>(                                \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04, hqvs); \
      error_estimate = result[(degree)];                                       \
      return result;                                                         \
    }                                                                        \
  }

// The error estimate must be computed in the Чебышёв basis because the elements
// of the monomial basis are not bounded.
#define PRINCIPIA_NEWHALL_MONOMIAL_APPROXIMATOR_SPECIALIZATION(degree)        \
  template<typename Value>                                                    \
  struct NewhallMonomialApproximator<Value, (degree)> {                       \
    static std::array<Value, (degree) + 1> HomogeneousCoefficients(           \
        HomogeneousQVs<Value> const& hqvs,                                    \
        Value& error_estimate) {                                              \
      error_estimate = DotProduct(                                            \
          newhall_c_matrix_чебышёв_degree_##degree##_divisions_8_w04          \
              .row<(degree)>(),                                               \
          hqvs);                                                              \
      return Multiply<(degree)>(                                              \
          newhall_c_matrix_monomial_degree_##degree##_divisions_8_w04, hqvs); \
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
  //TODO(phl)invert?
  HomogeneousQVs<Difference<Value>> hqvs;
  for (int i = 0, j = divisions; i < divisions + 1 && j >= 0; ++i, --j) {
    auto const& [q, v] = qvs[i];
    hqvs[j] = {q, v * duration_over_two};
  }

  auto const coefficients =
      NewhallЧебышёвApproximator<Difference<Value>, degree>::
          HomogeneousCoefficients(hqvs, error_estimate);
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
  HomogeneousQVs<Difference<Value>> hqvs;
  for (int i = 0, j = divisions;
       i < divisions + 1 && j >= 0;
       ++i, --j) {
    auto const& [q, v] = qvs[i];
    hqvs[j] = {q - origin, v * duration_over_two};
  }

  Instant const t_mid = Barycentre({t_min, t_max});
  return origin + Dehomogeneize<Difference<Value>, degree, Evaluator>(
                      NewhallMonomialApproximator<Difference<Value>, degree>::
                          HomogeneousCoefficients(hqvs, error_estimate),
                      /*scale=*/1.0 / duration_over_two,
                      t_mid);
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
      break;
  }
}

#undef PRINCIPIA_NEWHALL_APPROXIMATION_IN_MONOMIAL_BASIS_CASE

}  // namespace internal
}  // namespace _newhall
}  // namespace numerics
}  // namespace principia
