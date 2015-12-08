
#include "numerics/чебышёв_series.hpp"

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/serialization.hpp"
#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/newhall.mathematica.cpp"

namespace principia {

using geometry::DoubleOrQuantityOrMultivectorSerializer;
using geometry::Multivector;
using geometry::R3Element;

namespace numerics {

template<typename Vector>
ЧебышёвSeries<Vector>::ЧебышёвSeries(std::vector<Vector> const& coefficients,
                                     Instant const& t_min,
                                     Instant const& t_max)
    : coefficients_(coefficients),
      degree_(static_cast<int>(coefficients_.size()) - 1),
      t_min_(t_min),
      t_max_(t_max) {
  CHECK_LE(0, degree_) << "Degree must be at least 0";
  CHECK_LT(t_min_, t_max_) << "Time interval must not be empty";
  // Precomputed to save operations at the expense of some accuracy loss.
  Time const duration = t_max_ - t_min_;
  t_mean_ = t_min_ + 0.5 * duration;
  two_over_duration_ = 2 / duration;
}


template<typename Vector>
ЧебышёвSeries<Vector>::ЧебышёвSeries(ЧебышёвSeries&& other)
    : coefficients_(std::move(other.coefficients_)),
      degree_(std::move(other.degree_)),
      t_min_(std::move(other.t_min_)),
      t_max_(std::move(other.t_max_)),
      t_mean_(std::move(other.t_mean_)),
      two_over_duration_(std::move(other.two_over_duration_)) {}

template<typename Vector>
ЧебышёвSeries<Vector>& ЧебышёвSeries<Vector>::operator=(
    ЧебышёвSeries&& other) {
  coefficients_ = std::move(other.coefficients_);
  degree_ = other.degree_;
  t_min_ = std::move(other.t_min_);
  t_max_ = std::move(other.t_max_);
  t_mean_ = std::move(other.t_mean_);
  two_over_duration_ = std::move(other.two_over_duration_);
  return *this;
}

template<typename Vector>
bool ЧебышёвSeries<Vector>::operator==(ЧебышёвSeries const& right) const {
  return coefficients_ == right.coefficients_ &&
         t_min_ == right.t_min_ &&
         t_max_ == right.t_max_;
}

template<typename Vector>
bool ЧебышёвSeries<Vector>::operator!=(ЧебышёвSeries const& right) const {
  return !ЧебышёвSeries<Vector>::operator==(right);
}

template<typename Vector>
Instant const& ЧебышёвSeries<Vector>::t_min() const {
  return t_min_;
}

template<typename Vector>
Instant const& ЧебышёвSeries<Vector>::t_max() const {
  return t_max_;
}

template<typename Vector>
Vector const& ЧебышёвSeries<Vector>::last_coefficient() const {
  return coefficients_[degree_];
}

template<typename Vector>
class WTF {
 public:
  static Vector EvaluateImplementation(std::vector<Vector> const& coefficients,
                                       int const degree,
                                       double const scaled_t) {
    double const two_scaled_t = scaled_t + scaled_t;
    // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
    // computation errors.  But if it goes too far, something is broken.
    // TODO(phl): This should use DCHECK but these macros don't work because the
    // Principia projects don't define NDEBUG.
  #ifdef _DEBUG
    CHECK_LE(scaled_t, 1.1);
    CHECK_GE(scaled_t, -1.1);
  #endif

    Vector const c0 = coefficients[0];
    switch (degree) {
      case 0:
        return c0;
      case 1:
        return c0 + scaled_t * coefficients[1];
      default:
        Vector b_kplus2 = coefficients[degree];
        Vector b_kplus1 = coefficients[degree - 1] + two_scaled_t * b_kplus2;
        Vector b_k;
        for (int k = degree - 2; k >= 1; --k) {
          b_k = coefficients[k] + two_scaled_t * b_kplus1 - b_kplus2;
          b_kplus2 = b_kplus1;
          b_kplus1 = b_k;
        }
        return c0 + scaled_t * b_kplus1 - b_kplus2;
    }
  }
};

template<typename Scalar, typename Frame, int rank>
class WTF<Multivector<Scalar, Frame, rank>> {
 public:
  static Multivector<Scalar, Frame, rank> EvaluateImplementation(
      std::vector<Multivector<Scalar, Frame, rank>> const& coefficients,
      int const degree,
      double const scaled_t) {
    double const two_scaled_t = scaled_t + scaled_t;
    // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
    // computation errors.  But if it goes too far, something is broken.
    // TODO(phl): This should use DCHECK but these macros don't work because the
    // Principia projects don't define NDEBUG.
  #ifdef _DEBUG
    CHECK_LE(scaled_t, 1.1);
    CHECK_GE(scaled_t, -1.1);
  #endif

    R3Element<Scalar> const c0 = coefficients[0].coordinates();
    Scalar const c0x = c0.x;
    Scalar const c0y = c0.y;
    Scalar const c0z = c0.z;
    switch (degree) {
      case 0: {
        return Multivector<Scalar, Frame, rank>({c0x, c0y, c0z});
      }
      case 1: {
        R3Element<Scalar> const c1 = coefficients[1].coordinates();
        Scalar const c1x = c1.x;
        Scalar const c1y = c1.y;
        Scalar const c1z = c1.z;
        return Multivector<Scalar, Frame, rank>({c0x + scaled_t * c1x,
                                                 c0y + scaled_t * c1y,
                                                 c0z + scaled_t * c1z});
      }
      default: {
        R3Element<Scalar> const cd = coefficients[degree].coordinates();
        Scalar b_kplus2x = cd.x;
        Scalar b_kplus2y = cd.y;
        Scalar b_kplus2z = cd.z;
        R3Element<Scalar> const cdm1 = coefficients[degree - 1].coordinates();
        Scalar b_kplus1x = cdm1.x + two_scaled_t * b_kplus2x;
        Scalar b_kplus1y = cdm1.y + two_scaled_t * b_kplus2y;
        Scalar b_kplus1z = cdm1.z + two_scaled_t * b_kplus2z;
        Scalar b_k;
        for (int k = degree - 2; k >= 1; --k) {
          R3Element<Scalar> const ck = coefficients[k].coordinates();
          b_k = ck.x + two_scaled_t * b_kplus1x - b_kplus2x;
          b_kplus2x = b_kplus1x;
          b_kplus1x = b_k;
          b_k = ck.y + two_scaled_t * b_kplus1y - b_kplus2y;
          b_kplus2y = b_kplus1y;
          b_kplus1y = b_k;
          b_k = ck.z + two_scaled_t * b_kplus1z - b_kplus2z;
          b_kplus2z = b_kplus1z;
          b_kplus1z = b_k;
        }
        return Multivector<Scalar, Frame, rank>(
            {c0x + scaled_t * b_kplus1x - b_kplus2x,
             c0y + scaled_t * b_kplus1y - b_kplus2y,
             c0z + scaled_t * b_kplus1z - b_kplus2z});
      }
    }
  }
};

template<typename Vector>
Vector ЧебышёвSeries<Vector>::Evaluate(Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  return WTF<Vector>::EvaluateImplementation(coefficients_, degree_, scaled_t);
  //#define CLENSHAW_STEP1(n) \
//  bnplus2 = coefficients_[n];
//#define CLENSHAW_STEP2(n, nplus1) \
//  bnplus1 = coefficients_[n] + two_scaled_t * bnplus2;
//#define CLENSHAW_STEP3(n, nplus1, nplus2) \
//  bn = coefficients_[n] + two_scaled_t * bnplus1 - bnplus2;  \
//  bnplus2 = bnplus1; \
//  bnplus1 = bn;
}

template<typename Vector>
Variation<Vector> ЧебышёвSeries<Vector>::EvaluateDerivative(
    Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  double const two_scaled_t = scaled_t + scaled_t;
  // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  // TODO(phl): See above.
#ifdef _DEBUG
  CHECK_LE(scaled_t, 1.1);
  CHECK_GE(scaled_t, -1.1);
#endif

  Vector b_kplus2_vector{};
  Vector b_kplus1_vector{};
  Vector* b_kplus2 = &b_kplus2_vector;
  Vector* b_kplus1 = &b_kplus1_vector;
  Vector* const& b_k = b_kplus2;  // An overlay.
  for (int k = degree_ - 1; k >= 1; --k) {
    *b_k = coefficients_[k + 1] * (k + 1) +
           two_scaled_t * *b_kplus1 - *b_kplus2;
    Vector* const last_b_k = b_k;
    b_kplus2 = b_kplus1;
    b_kplus1 = last_b_k;
  }
  return (coefficients_[1] + two_scaled_t * *b_kplus1 - *b_kplus2) *
             two_over_duration_;
}

template<typename Vector>
void ЧебышёвSeries<Vector>::WriteToMessage(
    not_null<serialization::ЧебышёвSeries*> const message) const {
  using Serializer = DoubleOrQuantityOrMultivectorSerializer<
                          Vector,
                          serialization::ЧебышёвSeries::Coefficient>;

  for (auto const& coefficient : coefficients_) {
    Serializer::WriteToMessage(coefficient, message->add_coefficient());
  }
  t_min_.WriteToMessage(message->mutable_t_min());
  t_max_.WriteToMessage(message->mutable_t_max());
}

template<typename Vector>
ЧебышёвSeries<Vector> ЧебышёвSeries<Vector>::ReadFromMessage(
    serialization::ЧебышёвSeries const& message) {
  using Serializer = DoubleOrQuantityOrMultivectorSerializer<
                          Vector,
                          serialization::ЧебышёвSeries::Coefficient>;

  std::vector<Vector> coefficients;
  coefficients.reserve(message.coefficient_size());
  for (auto const& coefficient : message.coefficient()) {
    coefficients.push_back(Serializer::ReadFromMessage(coefficient));
  }
  return ЧебышёвSeries(coefficients,
                       Instant::ReadFromMessage(message.t_min()),
                       Instant::ReadFromMessage(message.t_max()));
}

template<typename Vector>
ЧебышёвSeries<Vector> ЧебышёвSeries<Vector>::NewhallApproximation(
    int const degree,
    std::vector<Vector> const& q,
    std::vector<Variation<Vector>> const& v,
    Instant const& t_min,
    Instant const& t_max) {
  // Only supports 8 divisions for now.
  int const kDivisions = 8;
  CHECK_EQ(kDivisions + 1, q.size());
  CHECK_EQ(kDivisions + 1, v.size());

  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Vector, 2 * kDivisions + 2> qv;
  for (int i = 0, j = 2 * kDivisions;
       i < kDivisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i];
    qv[j + 1] = v[i] * duration_over_two;
  }

  std::vector<Vector> coefficients;
  coefficients.reserve(degree);
  switch (degree) {
    case 3:
      coefficients = newhall_c_matrix_degree_3_divisions_8_w04 * qv;
      break;
    case 4:
      coefficients = newhall_c_matrix_degree_4_divisions_8_w04 * qv;
      break;
    case 5:
      coefficients = newhall_c_matrix_degree_5_divisions_8_w04 * qv;
      break;
    case 6:
      coefficients = newhall_c_matrix_degree_6_divisions_8_w04 * qv;
      break;
    case 7:
      coefficients = newhall_c_matrix_degree_7_divisions_8_w04 * qv;
      break;
    case 8:
      coefficients = newhall_c_matrix_degree_8_divisions_8_w04 * qv;
      break;
    case 9:
      coefficients = newhall_c_matrix_degree_9_divisions_8_w04 * qv;
      break;
    case 10:
      coefficients = newhall_c_matrix_degree_10_divisions_8_w04 * qv;
      break;
    case 11:
      coefficients = newhall_c_matrix_degree_11_divisions_8_w04 * qv;
      break;
    case 12:
      coefficients = newhall_c_matrix_degree_12_divisions_8_w04 * qv;
      break;
    case 13:
      coefficients = newhall_c_matrix_degree_13_divisions_8_w04 * qv;
      break;
    case 14:
      coefficients = newhall_c_matrix_degree_14_divisions_8_w04 * qv;
      break;
    case 15:
      coefficients = newhall_c_matrix_degree_15_divisions_8_w04 * qv;
      break;
    case 16:
      coefficients = newhall_c_matrix_degree_16_divisions_8_w04 * qv;
      break;
    case 17:
      coefficients = newhall_c_matrix_degree_17_divisions_8_w04 * qv;
      break;
    default:
      LOG(FATAL) << "Unexpected degree " << degree;
      break;
  }
  CHECK_EQ(degree + 1, coefficients.size());
  return ЧебышёвSeries(coefficients, t_min, t_max);
}

}  // namespace numerics
}  // namespace principia
