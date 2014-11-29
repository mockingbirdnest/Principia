#pragma once

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
inline R3Element<Scalar>::R3Element() : x(), y(), z() {}

template<typename Scalar>
inline R3Element<Scalar>::R3Element(Scalar const& x,
                                    Scalar const& y,
                                    Scalar const& z) : x(x), y(y), z(z) {}

template<typename Scalar>
inline Scalar& R3Element<Scalar>::operator[](int const index) {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      LOG(FATAL) << __FUNCSIG__ << ": index = " << index;
      noreturn();
  }
}

template<typename Scalar>
inline Scalar const& R3Element<Scalar>::operator[](int const index) const {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      LOG(FATAL) << __FUNCSIG__ << ": index = " << index;
      noreturn();
  }
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator+=(
    R3Element<Scalar> const& right) {
  return *this = *this + right;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator-=(
    R3Element<Scalar> const& right) {
  return *this = *this - right;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator*=(double const right) {
  return *this = *this * right;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator/=(double const right) {
  return *this = *this / right;
}

template<typename Scalar>
inline Scalar R3Element<Scalar>::Norm() const {
  return quantities::Sqrt(Dot(*this, *this));
}

template<typename Scalar>
inline R3Element<Scalar> operator+(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(+right.x, +right.y, +right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator-(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(-right.x, -right.y, -right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator+(
    R3Element<Scalar> const& left,
    R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator-(
    R3Element<Scalar> const& left,
    R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator*(double const left,
                                   R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left * right.x,
                           left * right.y,
                           left * right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                                   double const right) {
  return R3Element<Scalar>(left.x * right,
                           left.y * right,
                           left.z * right);
}

template<typename Scalar>
inline R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                                   double const right) {
  return R3Element<Scalar>(left.x / right,
                           left.y / right,
                           left.z / right);
}

template<typename LDimension, typename RScalar>
inline R3Element<quantities::Product<quantities::Quantity<LDimension>, RScalar>>
operator*(quantities::Quantity<LDimension> const& left,
          R3Element<RScalar> const& right) {
  return R3Element<quantities::Product<quantities::Quantity<LDimension>,
                                       RScalar>>(
      left * right.x,
      left * right.y,
      left * right.z);
}

template<typename LScalar, typename RDimension>
inline R3Element<quantities::Product<LScalar, quantities::Quantity<RDimension>>>
operator*(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3Element<quantities::Product<LScalar,
                                       quantities::Quantity<RDimension>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}

template<typename LScalar, typename RDimension>
inline R3Element<quantities::Quotient<LScalar,
                                      quantities::Quantity<RDimension>>>
operator/(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3Element<quantities::Quotient<LScalar,
                                        quantities::Quantity<RDimension>>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename Scalar>
bool operator==(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template<typename Scalar>
bool operator!=(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right) {
  return left.x != right.x || left.y != right.y || left.z != right.z;
}

template<typename Scalar>
std::string DebugString(R3Element<Scalar> const& r3_element) {
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
                         R3Element<Scalar> const& r3_element) {
  out << DebugString(r3_element);
  return out;
}

template<typename LScalar, typename RScalar>
inline R3Element<quantities::Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3Element<quantities::Product<LScalar, RScalar>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LScalar, typename RScalar>
inline quantities::Product<LScalar, RScalar> Dot(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace geometry
}  // namespace principia
