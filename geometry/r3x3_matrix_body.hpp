#pragma once

#include "geometry/r3x3_matrix.hpp"

#include <string>
#include <utility>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {

inline R3x3Matrix::R3x3Matrix(R3Element<double> const& row_x,
                              R3Element<double> const& row_y,
                              R3Element<double> const& row_z)
    : row_x_(row_x), row_y_(row_y), row_z_(row_z) {}

inline double R3x3Matrix::Trace() const {
  return row_x_.x + row_y_.y + row_z_.z;
}

FORCE_INLINE double R3x3Matrix::operator()(int const r, int const c) const {
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

inline R3x3Matrix R3x3Matrix::Transpose() const {
  return R3x3Matrix({row_x_.x, row_y_.x, row_z_.x},
                    {row_x_.y, row_y_.y, row_z_.y},
                    {row_x_.z, row_y_.z, row_z_.z});
}

inline R3x3Matrix& R3x3Matrix::operator+=(
    R3x3Matrix const& right) {
  return *this = *this + right;
}

inline R3x3Matrix& R3x3Matrix::operator-=(
    R3x3Matrix const& right) {
  return *this = *this - right;
}

inline R3x3Matrix& R3x3Matrix::operator*=(
    R3x3Matrix const& right) {
  return *this = *this * right;
}

inline R3x3Matrix& R3x3Matrix::operator*=(double const right) {
  return *this = *this * right;
}

inline R3x3Matrix& R3x3Matrix::operator/=(double const right) {
  return *this = *this / right;
}

inline R3x3Matrix R3x3Matrix::Identity() {
  return R3x3Matrix({1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1});
}

inline void R3x3Matrix::WriteToMessage(
    not_null<serialization::R3x3Matrix*> const message) const {
  row_x_.WriteToMessage(message->mutable_row_x());
  row_y_.WriteToMessage(message->mutable_row_y());
  row_z_.WriteToMessage(message->mutable_row_z());
}

inline R3x3Matrix R3x3Matrix::ReadFromMessage(
    serialization::R3x3Matrix const& message) {
  return R3x3Matrix(R3Element<double>::ReadFromMessage(message.row_x()),
                    R3Element<double>::ReadFromMessage(message.row_y()),
                    R3Element<double>::ReadFromMessage(message.row_z()));
}

inline R3x3Matrix operator+(R3x3Matrix const& right) {
  return R3x3Matrix(+right.row_x_, +right.row_y_, +right.row_z_);
}

inline R3x3Matrix operator-(R3x3Matrix const& right) {
  return R3x3Matrix(-right.row_x_, -right.row_y_, -right.row_z_);
}

inline R3x3Matrix operator+(R3x3Matrix const& left,
                            R3x3Matrix const& right) {
  return R3x3Matrix(left.row_x_ + right.row_x_,
                    left.row_y_ + right.row_y_,
                    left.row_z_ + right.row_z_);
}

inline R3x3Matrix operator-(R3x3Matrix const& left,
                            R3x3Matrix const& right) {
  return R3x3Matrix(left.row_x_ - right.row_x_,
                    left.row_y_ - right.row_y_,
                    left.row_z_ - right.row_z_);
}

inline R3x3Matrix operator*(R3x3Matrix const& left,
                            R3x3Matrix const& right) {
  R3x3Matrix const t_right = right.Transpose();
  return R3x3Matrix({Dot(left.row_x_, t_right.row_x_),
                     Dot(left.row_x_, t_right.row_y_),
                     Dot(left.row_x_, t_right.row_z_)},
                    {Dot(left.row_y_, t_right.row_x_),
                     Dot(left.row_y_, t_right.row_y_),
                     Dot(left.row_y_, t_right.row_z_)},
                    {Dot(left.row_z_, t_right.row_x_),
                     Dot(left.row_z_, t_right.row_y_),
                     Dot(left.row_z_, t_right.row_z_)});
}

template<typename Scalar>
R3Element<Scalar> operator*(R3x3Matrix const& left,
                            R3Element<Scalar> const& right) {
  return R3Element<Scalar>({Dot(left.row_x_, right),
                            Dot(left.row_y_, right),
                            Dot(left.row_z_, right)});
}

template<typename Scalar>
R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                            R3x3Matrix const& right) {
  R3x3Matrix const t_right = right.Transpose();
  return R3Element<Scalar>({Dot(left, t_right.row_x_),
                            Dot(left, t_right.row_y_),
                            Dot(left, t_right.row_z_)});
}


inline R3x3Matrix operator*(double const left,
                            R3x3Matrix const& right) {
  return R3x3Matrix(left * right.row_x_,
                    left * right.row_y_,
                    left * right.row_z_);
}

inline R3x3Matrix operator*(R3x3Matrix const& left,
                            double const right) {
  return R3x3Matrix(left.row_x_ * right,
                    left.row_y_ * right,
                    left.row_z_ * right);
}

inline R3x3Matrix operator/(R3x3Matrix const& left,
                            double const right) {
  return R3x3Matrix(left.row_x_ / right,
                    left.row_y_ / right,
                    left.row_z_ / right);
}

inline bool operator==(R3x3Matrix const& left,
                       R3x3Matrix const& right) {
  return left.row_x_ == right.row_x_ &&
         left.row_y_ == right.row_y_ &&
         left.row_z_ == right.row_z_;
}

inline bool operator!=(R3x3Matrix const& left,
                       R3x3Matrix const& right) {
  return left.row_x_ != right.row_x_ ||
         left.row_y_ != right.row_y_ ||
         left.row_z_ != right.row_z_;
}

inline std::string DebugString(R3x3Matrix const& r3x3_matrix) {
  std::string result = "{";
  result += DebugString(r3x3_matrix.row_x_);
  result += ", ";
  result += DebugString(r3x3_matrix.row_y_);
  result += ", ";
  result += DebugString(r3x3_matrix.row_z_);
  result +="}";
  return result;
}

inline std::ostream& operator<<(std::ostream& out,
                                R3x3Matrix const& r3x3_matrix) {
  out << DebugString(r3x3_matrix);
  return out;
}

}  // namespace geometry
}  // namespace principia
