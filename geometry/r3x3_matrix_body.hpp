
#pragma once

#include "geometry/r3x3_matrix.hpp"

#include <string>
#include <utility>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_r3x3_matrix {

using quantities::SIUnit;

template<typename Scalar>
R3x3Matrix<Scalar>::R3x3Matrix(R3Element<Scalar> const& row_x,
                               R3Element<Scalar> const& row_y,
                               R3Element<Scalar> const& row_z)
    : row_x_(row_x), row_y_(row_y), row_z_(row_z) {}

template<typename Scalar>
Scalar R3x3Matrix<Scalar>::Trace() const {
  return row_x_.x + row_y_.y + row_z_.z;
}

template<typename Scalar>
FORCE_INLINE(inline) Scalar R3x3Matrix<Scalar>::operator()(
    int const r, int const c) const {
  switch (r) {
    case 0:
      return row_x_[c];
    case 1:
      return row_y_[c];
    case 2:
      return row_z_[c];
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE
                  << " indices = {" << r << ", " << c << "}";
      base::noreturn();
  }
}

template<typename Scalar>
R3x3Matrix<Scalar> R3x3Matrix<Scalar>::Transpose() const {
  return R3x3Matrix({row_x_.x, row_y_.x, row_z_.x},
                    {row_x_.y, row_y_.y, row_z_.y},
                    {row_x_.z, row_y_.z, row_z_.z});
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator+=(
    R3x3Matrix const& right) {
  return *this = *this + right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator-=(
    R3x3Matrix const& right) {
  return *this = *this - right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator*=(
    R3x3Matrix const& right) {
  return *this = *this * right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator*=(double const right) {
  return *this = *this * right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator/=(double const right) {
  return *this = *this / right;
}

template<typename Scalar>
R3x3Matrix<Scalar> R3x3Matrix<Scalar>::Identity() {
  constexpr Scalar one = SIUnit<Scalar>();
  constexpr Scalar zero;
  return R3x3Matrix({one, zero, zero},
                    {zero, one, zero},
                    {zero, zero, one});
}

template<typename Scalar>
void R3x3Matrix<Scalar>::WriteToMessage(
    not_null<serialization::R3x3Matrix*> const message) const {
  row_x_.WriteToMessage(message->mutable_row_x());
  row_y_.WriteToMessage(message->mutable_row_y());
  row_z_.WriteToMessage(message->mutable_row_z());
}

template<typename Scalar>
R3x3Matrix<Scalar> R3x3Matrix<Scalar>::ReadFromMessage(
    serialization::R3x3Matrix const& message) {
  return R3x3Matrix(R3Element<double>::ReadFromMessage(message.row_x()),
                    R3Element<double>::ReadFromMessage(message.row_y()),
                    R3Element<double>::ReadFromMessage(message.row_z()));
}

template<typename Scalar>
R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix(+right.row_x_, +right.row_y_, +right.row_z_);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix(-right.row_x_, -right.row_y_, -right.row_z_);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& left,
                             R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix(left.row_x_ + right.row_x_,
                    left.row_y_ + right.row_y_,
                    left.row_z_ + right.row_z_);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& left,
                             R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix(left.row_x_ - right.row_x_,
                    left.row_y_ - right.row_y_,
                    left.row_z_ - right.row_z_);
}

template<typename LScalar, typename RScalar>
R3x3Matrix<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    R3x3Matrix<RScalar> const& right) {
  R3x3Matrix const t_right = right.Transpose();
  return R3x3Matrix<Product<LScalar, RScalar>>(
             {Dot(left.row_x_, t_right.row_x_),
              Dot(left.row_x_, t_right.row_y_),
              Dot(left.row_x_, t_right.row_z_)},
             {Dot(left.row_y_, t_right.row_x_),
              Dot(left.row_y_, t_right.row_y_),
              Dot(left.row_y_, t_right.row_z_)},
             {Dot(left.row_z_, t_right.row_x_),
              Dot(left.row_z_, t_right.row_y_),
              Dot(left.row_z_, t_right.row_z_)});
}

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3Element<Product<LScalar, RScalar>>({Dot(left.row_x_, right),
                                               Dot(left.row_y_, right),
                                               Dot(left.row_z_, right)});
}

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    R3Element<LScalar> const& left,
    R3x3Matrix<RScalar> const& right) {
  R3x3Matrix const t_right = right.Transpose();
  return R3Element<Product<LScalar, RScalar>>({Dot(left, t_right.row_x_),
                                               Dot(left, t_right.row_y_),
                                               Dot(left, t_right.row_z_)});
}


template<typename Scalar>
R3x3Matrix<Scalar> operator*(double const left,
                             R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix(left * right.row_x_,
                    left * right.row_y_,
                    left * right.row_z_);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator*(R3x3Matrix<Scalar> const& left,
                             double const right) {
  return R3x3Matrix(left.row_x_ * right,
                    left.row_y_ * right,
                    left.row_z_ * right);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator/(R3x3Matrix<Scalar> const& left,
                             double const right) {
  return R3x3Matrix(left.row_x_ / right,
                    left.row_y_ / right,
                    left.row_z_ / right);
}

template<typename LScalar, typename RScalar>
R3x3Matrix<Product<LScalar, RScalar>> KroneckerProduct(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3x3Matrix(left.x * right,
                    left.y * right,
                    left.z * right);
}

template<typename Scalar>
bool operator==(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right) {
  return left.row_x_ == right.row_x_ &&
         left.row_y_ == right.row_y_ &&
         left.row_z_ == right.row_z_;
}

template<typename Scalar>
bool operator!=(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right) {
  return left.row_x_ != right.row_x_ ||
         left.row_y_ != right.row_y_ ||
         left.row_z_ != right.row_z_;
}

template<typename Scalar>
std::string DebugString(R3x3Matrix<Scalar> const& r3x3_matrix) {
  std::string result = "{";
  result += DebugString(r3x3_matrix.row_x_);
  result += ", ";
  result += DebugString(r3x3_matrix.row_y_);
  result += ", ";
  result += DebugString(r3x3_matrix.row_z_);
  result +="}";
  return result;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix<Scalar> const& r3x3_matrix) {
  out << DebugString(r3x3_matrix);
  return out;
}

}  // namespace internal_r3x3_matrix
}  // namespace geometry
}  // namespace principia
