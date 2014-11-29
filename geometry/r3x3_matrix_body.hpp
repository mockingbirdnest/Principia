#pragma once

#include "geometry/r3x3_matrix.hpp"

#include <string>

#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {

inline R3x3Matrix::R3x3Matrix()
    : row_x_({1, 0, 0}),
      row_y_({0, 1, 0}),
      row_z_({0, 0, 1}) {}

inline R3x3Matrix::R3x3Matrix(R3Element<double> const& row_x,
                              R3Element<double> const& row_y,
                              R3Element<double> const& row_z)
    : row_x_(row_x), row_y_(row_y), row_z_(row_z) {}

R3x3Matrix R3x3Matrix::Transpose() const {
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

inline R3x3Matrix operator*(
    R3x3Matrix const& left,
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

bool operator==(R3x3Matrix const& left,
                R3x3Matrix const& right) {
  return left.row_x_ == right.row_x_ &&
         left.row_y_ == right.row_y_ &&
         left.row_z_ == right.row_z_;
}

bool operator!=(R3x3Matrix const& left,
                R3x3Matrix const& right) {
  return left.row_x_ != right.row_x_ ||
         left.row_y_ != right.row_y_ ||
         left.row_z_ != right.row_z_;
}

std::string DebugString(R3x3Matrix const& r3x3_matrix) {
  std::string result = "{";
  result += DebugString(r3x3_matrix.row_x_);
  result += ", ";
  result += DebugString(r3x3_matrix.row_y_);
  result += ", ";
  result += DebugString(r3x3_matrix.row_z_);
  result +="}";
  return result;
}

std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix const& r3x3_matrix) {
  out << DebugString(r3x3_matrix);
  return out;
}

}  // namespace geometry
}  // namespace principia
