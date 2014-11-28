#pragma once

#include "geometry/r3x3_matrix.hpp"

#include <assert.h>
#include <stdlib.h>

#include <string>

#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {

namespace {
// clang understands that this function is never used, but thinks that control
// reaches beyond |LOG(FATAL)| if it is not there.
// MSVC doesn't understand anything.
#ifdef __clang__
__attribute__((unused))
#endif
__declspec(noreturn) void noreturn() { exit(0); }
}  // namespace

inline R3x3Matrix::R3x3Matrix()
    : column_x_({1, 0, 0}),
      column_y_({0, 1, 0}),
      column_z_({0, 0, 1}) {}

inline R3x3Matrix::R3x3Matrix(R3Element<double> const& column_x,
                              R3Element<double> const& column_y,
                              R3Element<double> const& column_z)
    : column_x_(column_x), column_y_(column_y), column_z_(column_z) {}

inline R3x3Matrix& R3x3Matrix::operator+=(
    R3x3Matrix const& right) {
  return *this = *this + right;
}

inline R3x3Matrix& R3x3Matrix::operator-=(
    R3x3Matrix const& right) {
  return *this = *this - right;
}

inline R3x3Matrix& R3x3Matrix::operator*=(double const right) {
  return *this = *this * right;
}

inline R3x3Matrix& R3x3Matrix::operator/=(double const right) {
  return *this = *this / right;
}

inline R3x3Matrix operator+(R3x3Matrix const& right) {
  return R3x3Matrix(+right.column_x_, +right.column_y_, +right.column_z_);
}

inline R3x3Matrix operator-(R3x3Matrix const& right) {
  return R3x3Matrix(-right.column_x_, -right.column_y_, -right.column_z_);
}

inline R3x3Matrix operator+(
    R3x3Matrix const& left,
    R3x3Matrix const& right) {
  return R3x3Matrix(left.column_x_ + right.column_x_,
                           left.column_y_ + right.column_y_,
                           left.column_z_ + right.column_z_);
}

inline R3x3Matrix operator-(
    R3x3Matrix const& left,
    R3x3Matrix const& right) {
  return R3x3Matrix(left.column_x_ - right.column_x_,
                           left.column_y_ - right.column_y_,
                           left.column_z_ - right.column_z_);
}

inline R3x3Matrix operator*(double const left,
                                   R3x3Matrix const& right) {
  return R3x3Matrix(left * right.column_x_,
                           left * right.column_y_,
                           left * right.column_z_);
}

inline R3x3Matrix operator*(R3x3Matrix const& left,
                                   double const right) {
  return R3x3Matrix(left.column_x_ * right,
                           left.column_y_ * right,
                           left.column_z_ * right);
}

inline R3x3Matrix operator/(R3x3Matrix const& left,
                                   double const right) {
  return R3x3Matrix(left.column_x_ / right,
                           left.column_y_ / right,
                           left.column_z_ / right);
}

bool operator==(R3x3Matrix const& left,
                R3x3Matrix const& right) {
  return left.column_x_ == right.column_x_ &&
         left.column_y_ == right.column_y_ &&
         left.column_z_ == right.column_z_;
}

bool operator!=(R3x3Matrix const& left,
                R3x3Matrix const& right) {
  return left.column_x_ != right.column_x_ ||
         left.column_y_ != right.column_y_ ||
         left.column_z_ != right.column_z_;
}

std::string DebugString(R3x3Matrix const& r3_element) {
  std::string result = "{";
  result += DebugString(r3_element.column_x_);
  result += ", ";
  result += DebugString(r3_element.column_y_);
  result += ", ";
  result += DebugString(r3_element.column_z_);
  result +="}";
  return result;
}

std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix const& r3_element) {
  out << DebugString(r3_element);
  return out;
}

}  // namespace geometry
}  // namespace principia
