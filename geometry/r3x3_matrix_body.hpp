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

// We want zero initialization here, so the default constructor won't do.
inline R3x3Matrix::R3x3Matrix() : x(), y(), z() {}

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
  return R3x3Matrix(+right.x, +right.y, +right.z);
}

inline R3x3Matrix operator-(R3x3Matrix const& right) {
  return R3x3Matrix(-right.x, -right.y, -right.z);
}

inline R3x3Matrix operator+(
    R3x3Matrix const& left,
    R3x3Matrix const& right) {
  return R3x3Matrix(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
}

inline R3x3Matrix operator-(
    R3x3Matrix const& left,
    R3x3Matrix const& right) {
  return R3x3Matrix(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
}

inline R3x3Matrix operator*(double const left,
                                   R3x3Matrix const& right) {
  return R3x3Matrix(left * right.x,
                           left * right.y,
                           left * right.z);
}

inline R3x3Matrix operator*(R3x3Matrix const& left,
                                   double const right) {
  return R3x3Matrix(left.x * right,
                           left.y * right,
                           left.z * right);
}

inline R3x3Matrix operator/(R3x3Matrix const& left,
                                   double const right) {
  return R3x3Matrix(left.x / right,
                           left.y / right,
                           left.z / right);
}

bool operator==(R3x3Matrix const& left,
                R3x3Matrix const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool operator!=(R3x3Matrix const& left,
                R3x3Matrix const& right) {
  return left.x != right.x || left.y != right.y || left.z != right.z;
}

std::string DebugString(R3x3Matrix const& r3_element) {
  std::string result = "{";
  result += quantities::DebugString(r3_element.x);
  result += ", ";
  result += quantities::DebugString(r3_element.y);
  result += ", ";
  result += quantities::DebugString(r3_element.z);
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
