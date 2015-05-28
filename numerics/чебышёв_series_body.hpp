
#include "numerics/чебышёв_series.hpp"

#include "glog/logging.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/newhall.mathematica.cpp"
#include "quantities/serialization.hpp"

namespace principia {

using quantities::QuantityOrDoubleSerializer;

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
      degree_(other.degree_),
      t_min_(std::move(other.t_min_)),
      t_max_(std::move(other.t_max_)),
      t_mean_(std::move(other.t_mean_)),
      two_over_duration_(std::move(other.two_over_duration_)) {}

template<typename Vector>
ЧебышёвSeries<Vector>& ЧебышёвSeries<Vector>::operator=(ЧебышёвSeries&& other) {
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
Vector ЧебышёвSeries<Vector>::Evaluate(Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  double const two_scaled_t = scaled_t + scaled_t;
  // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  CHECK_LE(scaled_t, 1.1);
  CHECK_GE(scaled_t, -1.1);

  Vector b_kplus2{};
  Vector b_kplus1{};
  Vector b_k{};
  for (int k = degree_; k >= 1; --k) {
    b_k = coefficients_[k] + two_scaled_t * b_kplus1 - b_kplus2;
    b_kplus2 = b_kplus1;
    b_kplus1 = b_k;
  }
  return coefficients_[0] + scaled_t * b_kplus1 - b_kplus2;
}

template<typename Vector>
Variation<Vector> ЧебышёвSeries<Vector>::EvaluateDerivative(
    Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  double const two_scaled_t = scaled_t + scaled_t;
  // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  CHECK_LE(scaled_t, 1.1);
  CHECK_GE(scaled_t, -1.1);

  Vector b_kplus2{};
  Vector b_kplus1{};
  Vector b_k{};
  for (int k = degree_ - 1; k >= 1; --k) {
    b_k = coefficients_[k + 1] * (k + 1) + two_scaled_t * b_kplus1 - b_kplus2;
    b_kplus2 = b_kplus1;
    b_kplus1 = b_k;
  }
  return (coefficients_[1] + two_scaled_t * b_kplus1 - b_kplus2) *
             two_over_duration_;
}

template<typename Vector>
void ЧебышёвSeries<Vector>::WriteToMessage(
    not_null<serialization::ЧебышёвSeries*> const message) const {
  using Serializer =
      QuantityOrDoubleSerializer<Vector,
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
  using Serializer =
      QuantityOrDoubleSerializer<Vector,
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
