#pragma once

#include "numerics/newhall.hpp"

#include <vector>

#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_newhall {

using quantities::Exponentiation;
using quantities::Frequency;
using quantities::Pow;
using quantities::Time;

// Only supports 8 divisions for now.
constexpr int divisions = 8;

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Time, degree, Evaluator> Dehomogeneize(
    FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale) {
  using P = PolynomialInMonomialBasis<Vector, Time, degree, Evaluator>;
  typename P::Coefficients dehomogeneized_coefficients;
  Dehomogeneizer<Vector, typename P::Coefficients, degree, /*k=*/0>::Convert(
      homogeneous_coefficients,
      scale,
      dehomogeneized_coefficients);
  return P(dehomogeneized_coefficients);
}

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree, int k>
struct Dehomogeneizer {
  static void Convert(
      FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
      Frequency const& scale,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree>
struct Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, degree> {
  static void Convert(
      FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
      Frequency const& time_scale,
      DehomogeneizedCoefficients& dehomogeneized_coefficients);
};

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree, int k>
void Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, k>::Convert(
    FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
    Frequency const& scale,
    DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<k>(dehomogeneized_coefficients) =
      homogeneous_coefficients[k] * Pow<k>(scale);
  Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, k + 1>::Convert(
      homogeneous_coefficients,
      scale,
      dehomogeneized_coefficients);
}

template<typename Vector,
         typename DehomogeneizedCoefficients, int degree>
void Dehomogeneizer<Vector, DehomogeneizedCoefficients, degree, degree>::
Convert(FixedVector<Vector, degree + 1> const& homogeneous_coefficients,
        Frequency const& scale,
        DehomogeneizedCoefficients& dehomogeneized_coefficients) {
  std::get<degree>(dehomogeneized_coefficients) =
      homogeneous_coefficients[degree] * Pow<degree>(scale);
}

template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator {
  static FixedVector<Vector, degree + 1> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv);
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 3, Evaluator> {
  static FixedVector<Vector, 4> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_3_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 4, Evaluator> {
  static FixedVector<Vector, 5> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_4_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 5, Evaluator> {
  static FixedVector<Vector, 6> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_5_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 6, Evaluator> {
  static FixedVector<Vector, 7> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_6_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 7, Evaluator> {
  static FixedVector<Vector, 8> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_7_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 8, Evaluator> {
  static FixedVector<Vector, 9> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_8_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 9, Evaluator> {
  static FixedVector<Vector, 10> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_9_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 10, Evaluator> {
  static FixedVector<Vector, 11> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_10_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 11, Evaluator> {
  static FixedVector<Vector, 12> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_11_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 12, Evaluator> {
  static FixedVector<Vector, 13> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_12_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 13, Evaluator> {
  static FixedVector<Vector, 14> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_13_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 14, Evaluator> {
  static FixedVector<Vector, 15> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_14_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 15, Evaluator> {
  static FixedVector<Vector, 16> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_15_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 16, Evaluator> {
  static FixedVector<Vector, 17> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_16_divisions_8_w04 * qv;
  }
};

template<typename Vector,
         template<typename, typename, int> class Evaluator>
struct NewhallAppromixator<Vector, 17, Evaluator> {
  static FixedVector<Vector, 18> HomogeneousCoefficients(
      FixedVector<Vector, 2 * divisions + 2> const& qv) {
    return newhall_c_matrix_monomial_degree_17_divisions_8_w04 * qv;
  }
};

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

  // TODO(phl): Implement error estimation.
  error_estimate = Vector{};
  return Dehomogeneize<Vector, degree, Evaluator>(
             NewhallAppromixator<Vector, degree, Evaluator>::
                 HomogeneousCoefficients(qv),
             /*scale=*/1.0 / duration_over_two);
}

}  // namespace internal_newhall
}  // namespace numerics
}  // namespace principia
