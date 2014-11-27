#pragma once

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
template<typename Scalar>
inline R3x3Matrix<Scalar>::R3x3Matrix() : x(), y(), z() {}

template<typename Scalar>
inline R3x3Matrix<Scalar>::R3x3Matrix(R3Element<Scalar> const& column_x,
                                      R3Element<Scalar> const& column_y,
                                      R3Element<Scalar> const& column_z)
    : column_x_(column_x), column_y_(column_y), column_z_(column_z) {}

template<typename Scalar>
inline R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator+=(
    R3x3Matrix<Scalar> const& right) {
  return *this = *this + right;
}

template<typename Scalar>
inline R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator-=(
    R3x3Matrix<Scalar> const& right) {
  return *this = *this - right;
}

template<typename Scalar>
inline R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator*=(double const right) {
  return *this = *this * right;
}

template<typename Scalar>
inline R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator/=(double const right) {
  return *this = *this / right;
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix<Scalar>(+right.x, +right.y, +right.z);
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix<Scalar>(-right.x, -right.y, -right.z);
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator+(
    R3x3Matrix<Scalar> const& left,
    R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix<Scalar>(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator-(
    R3x3Matrix<Scalar> const& left,
    R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix<Scalar>(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator*(double const left,
                                   R3x3Matrix<Scalar> const& right) {
  return R3x3Matrix<Scalar>(left * right.x,
                           left * right.y,
                           left * right.z);
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator*(R3x3Matrix<Scalar> const& left,
                                   double const right) {
  return R3x3Matrix<Scalar>(left.x * right,
                           left.y * right,
                           left.z * right);
}

template<typename Scalar>
inline R3x3Matrix<Scalar> operator/(R3x3Matrix<Scalar> const& left,
                                   double const right) {
  return R3x3Matrix<Scalar>(left.x / right,
                           left.y / right,
                           left.z / right);
}

template<typename LDimension, typename RScalar>
inline R3x3Matrix<quantities::Product<quantities::Quantity<LDimension>, RScalar>>
operator*(quantities::Quantity<LDimension> const& left,
          R3x3Matrix<RScalar> const& right) {
  return R3x3Matrix<quantities::Product<quantities::Quantity<LDimension>,
                                       RScalar>>(
      left * right.x,
      left * right.y,
      left * right.z);
}

template<typename LScalar, typename RDimension>
inline R3x3Matrix<quantities::Product<LScalar, quantities::Quantity<RDimension>>>
operator*(R3x3Matrix<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3x3Matrix<quantities::Product<LScalar,
                                       quantities::Quantity<RDimension>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}

template<typename LScalar, typename RDimension>
inline R3x3Matrix<quantities::Quotient<LScalar,
                                      quantities::Quantity<RDimension>>>
operator/(R3x3Matrix<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3x3Matrix<quantities::Quotient<LScalar,
                                        quantities::Quantity<RDimension>>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename Scalar>
bool operator==(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template<typename Scalar>
bool operator!=(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right) {
  return left.x != right.x || left.y != right.y || left.z != right.z;
}

template<typename Scalar>
std::string DebugString(R3x3Matrix<Scalar> const& r3_element) {
  std::string result = "{";
  result += quantities::DebugString(r3_element.x);
  result += ", ";
  result += quantities::DebugString(r3_element.y);
  result += ", ";
  result += quantities::DebugString(r3_element.z);
  result +="}";
  return result;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix<Scalar> const& r3_element) {
  out << DebugString(r3_element);
  return out;
}

template<typename LScalar, typename RScalar>
inline R3x3Matrix<quantities::Product<LScalar, RScalar>> Cross(
    R3x3Matrix<LScalar> const& left,
    R3x3Matrix<RScalar> const& right) {
  return R3x3Matrix<quantities::Product<LScalar, RScalar>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LScalar, typename RScalar>
inline quantities::Product<LScalar, RScalar> Dot(
    R3x3Matrix<LScalar> const& left,
    R3x3Matrix<RScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace geometry
}  // namespace principia
